# A script to process ambient vibration records
from __future__ import print_function
# Use the print function to be able to switch easily between stdout and a file
from mpi4py import MPI
import os
import sys
import shutil
import time

from math import ceil
from obspy import read, Stream,  Trace, UTCDateTime
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import cheb2ord, cheby2, zpk2sos


from ants_2.tools import mergetraces as mt
from ants_2.tools import event_excluder as ee
from ants_2.tools.preprocess import Preprocess

from ants_2.config import ConfigPreprocess
cfg = ConfigPreprocess()

class Preprocessing(object):
    #ToDo find sensible names and task subdivisions
    

    def __init__(self):
        comm = MPI.COMM_WORLD
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()
        print("Hello from rank %g" %self.rank)
        print("Size is %g" %self.size)


    def preprocess(self):
        """
        
        This script preprocesses the MSEED files in the input directories 
        specified in the input file.
     
        
        """
    
        outdir = os.path.join('data','processed')
         
        try:
            os.mkdir(outdir)
        except OSError:
            pass
        
           
    #- check what input is, list input from different directories 
        
        content = self.find_rawdatafiles()
           
    # output file
        output_file = os.path.join(outdir,'processing_report_rank%g.txt' %self.rank)
        if os.path.exists(output_file):
            ofid = open(output_file,'a')
        else:
            ofid = open(output_file,'w')
            print('UPDATING, Date:',file=ofid)
            print(time.strftime('%Y.%m.%dT%H:%M'),file=ofid)
        #===============================================================================
        
        # Input file loop
        #===============================================================================
        rankdir = os.path.join(outdir,'rank_%g' %self.rank)
        
        if os.path.exists(rankdir)==False:
            os.mkdir(rankdir)
            
    
        
        for filepath in content[self.rank::self.size]:
            
            
            data = self.read_file(filepath,ofid)
            Fs_old = data[0].stats.sampling_rate
            id = data[0].id
            print(id)
            if data == []:
                continue
            
            
            #- clean the data merging segments with gap shorter than a specified number of seconds:
            data = mt.mergetraces(data,cfg.Fs_old,cfg.quality_maxgapsec,ofid)
            
            # splitting to avoid including gaps ToDo: Check if it is at all necessary
            data.split()
            
            #- initialize stream to 'recollect' the split traces
            colloc_data=Stream()
            
          
            #- slice traces into shorter segments
            data, ndata = self.slice_traces(data,ofid)
            
            #- trim to next second before decimation, to keep as synchronous as possible #ToDo: This should go into preprocess
            if cfg.wins_trim:
                data=self.trim_next_sec(data,ofid)
            
            # ToDo: Prettier event excluder
            if cfg.event_exclude:
                print(int(cfg.event_exclude))
                print('excluding events....')
                ee.event_exclude(data,windows,n_compare,min_freq,\
                factor_enrg=1.,taper_perc=0.05,thresh_stdv=1.,undo_taper=True)
            
            prep = Preprocess(data,cfg.verbose,ofid)
            
        
            
            prep.check_nan_inf()
            if cfg.wins_detrend:
                prep.detrend()
            if cfg.wins_demean:
                prep.demean()
            if cfg.wins_taper is not None:
                prep.taper(cfg.wins_taper)
            # ToDo cap glitches
            if cfg.Fs_new[-1] != Fs_old:
                # add antialias filter
                # Get filter coeff
                z, p, k = self.cheby2_lowpass(Fs_old,cfg.Fs_antialias_factor *
                 cfg.Fs_new[-1])
                sos = zpk2sos(z, p, k)
                prep.downsampling(sos)
                
            if cfg.instr_correction:
                # add inventory
                
                if cfg.instr_correction_input == 'resp':
                    respfile = os.path.join('meta','RESP.{}'.format(data[0].id))
                    inv = {
                        'filename': respfile, 
                        'units': cfg.instr_correction_unit
                                    }
                
                elif cfg.instr_correction_input == 'staxml':
                    sxmlfile = os.path.join('meta','{}.{}.xml'.format(data[0].stats.network,
                    data[0].stats.station))
                    inv = read_inventory(sxmlfile)
                                    
                prep.remove_response(inv,pre_filt=cfg.instr_correction_prefilt,
                waterlevel=cfg.instr_correction_waterlevel,output=cfg.instr_correction_unit)
            prep.check_nan_inf()
            
            
            data = prep.stream
            
            if len(data) == 0: 
                print('*** NO data returned from this file: '+filepath.split('/')[-1])
                continue
                
            data=mt.mergetraces(data,cfg.Fs_old,cfg.quality_maxgapsec,ofid)
            data._cleanup()
    
            for trace in data:
                
                
                filepathnew = getfilepath(rankdir,trace.stats)
                trace.write(filepathnew,format=trace.stats._format)
                           
                if cfg.verbose:
                    print('* renamed file: '+filepathnew,file=ofid)
            
            del data
            
        
        print("Rank %g has completed processing." %self.rank,file=ofid)
        ofid.close()
        
        try:
            os.system('mv '+rankdir+'/* '+rankdir+'/../')
        except:
            pass
        os.system('rmdir '+rankdir)
                
    def getfilepath(self,rankdir,stats,startonly=False):
        
        inf = [
            stats.network,
        stats.station,
        stats.location,
        stats.channel,
        stats._format
        ]
        
            
        t1=stats.starttime.strftime('%Y.%j.%H.%M.%S')
        t2=stats.endtime.strftime('%Y.%j.%H.%M.%S')
        if startonly: 
            t2 = '*'
        
        inf.append(t1,t2)
        
        #yr=str(t1[0:4])
        
        filenew = '{}.{}.{}.{}.{}.{}.{}'.format(inf)
        filepathnew = os.path.join(rankdir,filenew)
        #filepathnew=rankdir+'/'+network+'.'+station+'.'+location+'.'+\
        #channel+'.' + t1 + '.' +t2+'.'+prepname+'.'+format 
        
        return filepathnew
            
    
    def find_rawdatafiles(self):
        
        indirs = cfg.input_dirs
        format = cfg.input_format.lower()
        content=list()
        for indir in indirs:
            print(indir)
            content.extend(glob(os.path.join(indir,'*'+format)))
            content.extend(glob(os.path.join(indir,'*'+format.upper())))
            
        content.sort()
        print(content)
        return content
        
    def read_file(self,filepath,ofid):
        
        if cfg.verbose:
            print('===========================================================',\
            file=ofid)
            print('* opening file: '+filepath+'\n',file=ofid)
            
        #- read data
        try:
            data=read(filepath)
        except (TypeError, IOError):
            if cfg.verbose==True: print('** file wrong type or not found, skip.',
            file=ofid)
            return []
        except:
            if cfg.verbose: print('** unexpected read error, skip.',
            file=ofid)
            return []
            
        #- check if this file contains data
        if len(data) == 0:
            print('** file contains no data, skip.',file=ofid)
            return []
        
        return data
    
    
    
    def slice_traces(self,s, ofid):
        """
        Slice an ObsPy stream object with multiple traces; The stream of new 
        (sliced) traces merely contains references to the original trace.
        """
        length_in_sec = cfg.wins_len_sec
        min_len = cfg.quality_minlengthsec
        
        s_new=Stream()
        
        #- set initial start time
        starttime=s[0].stats.starttime
        
        
        for k in np.arange(len(s)):
               
            #- march through the trace until the endtime is reached
            while starttime < s[k].stats.endtime-min_len:
                
                s_part = s[k].slice(starttime,starttime+length_in_sec-1./
                (s[k].stats.sampling_rate))
                
                
                starttime += length_in_sec
                
        n_traces=len(s_new)
        
        if cfg.verbose:
            print('* contains %g trace(s)' %n_traces,file=ofid)
            
        return s_new, n_traces
        
    #ToDo: look into this
    def trim_next_sec(self,data,ofid):
        
        """ 
        Trim data to the full second. Ensures that recordings start and end with a full second.
    
        data=trim_next_sec(data,verbose)
    
        data:       Is an obspy stream or trace. The returned stream/trace is a bit shorter.
        verbose:    Talk or not.
        ofid: Output file in 'check' mode; otherwise, 'None' is the default inherited from main.
        
        """
    
        
    
        if isinstance(data,Trace):
            sec_to_remove=data.stats.starttime.microsecond/1e6
            sec_to_add=1.-sec_to_remove
            if sec_to_add > data.stats.delta:
                if cfg.verbose: 
                        print('* Trimming to full second.\n',file=ofid)
                data.trim(starttime=data.stats.starttime+sec_to_add,nearest_sample=True)
            
        elif isinstance(data,Stream):
            for tr in data:
                starttime=tr.stats.starttime
                sec_to_remove=tr.stats.starttime.microsecond/1e6
                sec_to_add=1.-sec_to_remove
                if sec_to_add > tr.stats.delta:
                    if cfg.verbose: 
                            print('* Trimming to full second.\n',file=ofid)
                    tr.trim(starttime=tr.stats.starttime+sec_to_add,nearest_sample=True)
                
        return data
        
    def cheby2_lowpass(self,df,freq,maxorder=8):
        # From obspy
        nyquist = df * 0.5
        # rp - maximum ripple of passband, rs - attenuation of stopband
        rp, rs, order = 1, 96, 1e99
        ws = freq / nyquist  # stop band frequency
        wp = ws  # pass band frequency
        # raise for some bad scenarios
        if ws > 1:
            ws = 1.0
            msg = "Selected corner frequency is above Nyquist. " + \
                  "Setting Nyquist as high corner."
            warnings.warn(msg)
        while True:
            if order <= maxorder:
                break
            wp = wp * 0.99
            order, wn = cheb2ord(wp, ws, rp, rs, analog=0)
        return cheby2(order, rs, wn, btype='low', analog=0, output='zpk')
        
    
 #   if check==True:
 #       ctr=trace.copy()
 #       ctr.stats.network='Original Data'
 #       ctr.stats.station=''
 #       ctr.stats.location=''
 #       ctr.stats.channel=''
 #       cstr=Stream(ctr)
 #       print(trace,file=dfile)
 #       dfile.write('-----------------------------------------------\n')
 #       dfile.write('Original\n')
 #       print(trace.data[0:20],file=dfile)
 #       dfile.write('\n')
 
#if check==True:
#    ctr = newtrace.copy()
#    ctr.stats.network='After IC to '+unit
#    ctr.stats.station=''
#    ctr.stats.location=''
#    ctr.stats.channel=''
#    
#    cstr.append(ctr)
#    cstr.plot(outfile=datadir+'/processed/out/'+\
#        filepath.split('/')[-1]+'.'+prepname+'.png',equal_scale=False)
#    cstr.trim(endtime=cstr[0].stats.starttime+3600)
#    cstr.plot(outfile=datadir+'/processed/out/'+\
#        filepath.split('/')[-1]+'.'+prepname+'.1hr.png',equal_scale=False)
#    dfile.write('Instrument response removed\n')
#    print(newtrace.data[0:20],file=dfile)
#    dfile.write('\n')
        #===============================================================================
        
        # trace loop
        #===============================================================================
#        
#        for itrace in range(ndata):
#            
#            # ToDo: put test run back in!
#            
#            trace = data[itrace]
#            
#            if trace.stats.npts / cfg.Fs_new[-1] < cfg.quality_minlengthsec:
#                continue
#            
#                
#            
#            
#            #if update == True:
#            #    if len(glob(getfilepath(rankdir,trace.stats,prepname,True))) > 0:
#            #        print('File already processed, proceeding...',file=ofid)
#            #        print(trace)
#            #        print('File already processed, proceeding...',file=None)
#            #        
#            #        break
#            #    else:
#            #        print('Updating...',file=ofid)
#            
#            if cfg.verbose: print('-----------------------------------------',
#            file=ofid)
#
#        
#            
#                
#            #- merge all into final trace   =========================================================
#            colloc_data+=newtrace
             
            #- flush buffer of output file ========================================================
#            ofid.flush()
#            
#            del newtrace
    