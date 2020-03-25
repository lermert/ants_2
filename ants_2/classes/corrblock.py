# Correlation block object:
from __future__ import print_function
from obspy import Stream, UTCDateTime, read
from scipy.signal import sosfilt
import numpy as np
import os
import re


from ants_2.tools.bookkeep import name_correlation_file
from ants_2.tools.util import get_geoinf
from ants_2.classes.corrtrace import CorrTrace
from ants_2.tools.correlations import cross_covar
from ants_2.tools.treatment import ram_norm, whiten, cap, bandpass
# list of possible channels combinations indicating that the data needs to be rotated.
horizontals = ['RR','RT','TR','TT','TZ','ZT','RZ','ZR']
import matplotlib.pyplot as plt


class CorrBlock(object):


# - initialize with station pairs
    def __init__(self,block,cfg):


        self.cfg = cfg
        self._correlations = {}


        self.inv = block.inventory
        self.channels = block.channels
        self.station_pairs = block.station_pairs
        self.channel_pairs = block.channel_pairs


        self.initialize_data()
        self.sampling_rate = self.data[0].stats.sampling_rate
        self.delta = self.data[0].stats.delta

        

        for cp in block.channel_pairs:
            for pair in cp:

                pair = [re.sub('E$','T',str) for str in pair]
                pair = [re.sub('N$','R',str) for str in pair]
                cp_name = '{}--{}'.format(*pair)
                preprstring = self.get_prepstring()

                
                self._correlations[cp_name] = CorrTrace(pair[0],pair[1],
                self.sampling_rate,stck_int=cfg.interm_stack,
                prepstring=self.get_prepstring(),
                window_length=cfg.time_window_length,
                overlap=cfg.time_overlap,corr_params=None)
                


        
        if any([i in cfg.corr_tensorcomponents for i in horizontals]):

            self.baz1 = []
            self.baz2 = []

            for pair in self.station_pairs:
                try:
                    geoinf = get_geoinf(pair[0],pair[1])
                    # - find azimuth, backazimuth
                    self.baz2.append(geoinf[5])
                    self.baz1.append(geoinf[6])
                except:
                    self.baz2.append(0)
                    self.baz1.append(0)

    

    def run(self,output_file=None):


        print('Working on station pairs:')
        for sta in self.station_pairs:
            print("{}--{}".format(sta[0],sta[1]))

        t_0 = UTCDateTime(self.cfg.time_begin)
        t_end = UTCDateTime(self.cfg.time_end)
        win_len_seconds = self.cfg.time_window_length
        win_len_samples = int(round(win_len_seconds*self.sampling_rate))
        min_len_samples = int(round(self.cfg.time_min_window*self.sampling_rate))
        max_lag_samples = int(round(self.cfg.corr_maxlag * self.sampling_rate))
        
        
        if self.cfg.bandpass is not None:
            fmin = self.cfg.bandpass[0]
            fmax = self.cfg.bandpass[1]
            if fmax <= fmin:
                msg = "Bandpass upper corner frequency must be above lower corner frequency."
                raise ValueError(msg)

            order = self.cfg.bandpass[2]
            sos = bandpass(freqmin=fmin,freqmax=fmax,
                df=self.sampling_rate,corners=order)


        # Time loop
        t = t_0
        print(t)    
        while t <= t_end - (win_len_seconds - self.delta):
            

            print(t,file=output_file)
            
            
            
            # - check endtime, if necessary, add data from 'later' file
            self.update_data(t, win_len_seconds)
            
            

            # - slice the traces
            
            if self.cfg.time_overlap == 0:
                windows = self.data.slice(t, t + win_len_seconds - self.delta)
            else:
                # - deepcopy is used so that processing is not applied directly on the data stream
                # - This is much more expensive than using non-overlapping windows, due to the copying
                windows = self.data.slice(t, t + win_len_seconds - self.delta).copy()
            
            

            
            # - Apply preprocessing
            for w in windows:
                
            
                if self.cfg.bandpass is not None:
                    w_temp = sosfilt(sos,w.data)
                    w.data = sosfilt(sos,w_temp[::-1])[::-1]

                if self.cfg.cap_glitch:
                    cap(w,self.cfg.cap_thresh)

                if self.cfg.whiten:
                    whiten(w,self.cfg.white_freqmin,
                        self.cfg.white_freqmax,
                        self.cfg.white_taper_samples)

                if self.cfg.onebit:
                    w.data = np.sign(w.data)

                if self.cfg.ram_norm:
                    ram_norm(w,self.cfg.ram_window,self.cfg.ram_prefilt)
                

            # - station pair loop
            for sp_i in range(len(self.station_pairs)):

                pair = self.station_pairs[sp_i]
                
                # - select traces
                [net1, sta1] = pair[0].split('.')
                [net2, sta2] = pair[1].split('.')
                
                str1 = windows.select(network=net1, station=sta1)
                str2 = windows.select(network=net2, station=sta2)
                

                # - if horizontal components are involved, copy and rotate
                if any([i in self.cfg.corr_tensorcomponents for i in horizontals]):

                    str1, str2 = self.rotate(str1,str2,self.baz1[sp_i],self.baz2[sp_i])
                    
                # - channel loop
                
                for cpair in self.channel_pairs[sp_i]:

                    cpair = [re.sub('E$','T',str) for str in cpair]
                    cpair = [re.sub('N$','R',str) for str in cpair]
                
                    cp_name = '{}--{}'.format(*cpair)
                    print(cp_name,file=output_file)
                    

                    loc1, cha1 = cpair[0].split('.')[2:4]
                    loc2, cha2 = cpair[1].split('.')[2:4]


                    try:
                        tr1 = str1.select(location=loc1,channel=cha1)[0]
                        tr2 = str2.select(location=loc2,channel=cha2)[0]
                        
                    except IndexError:
                        print("Channel not found",file=output_file)
                        continue
                        
                    # - check minimum length requirement
                    # - Quite often not fulfilled due to data gaps
                    if tr1.stats.npts < min_len_samples:
                        print("Trace length < min samples\n",file=output_file)
                        continue
                        
                    if tr2.stats.npts < min_len_samples:
                        print("Trace length < min samples\n",file=output_file)
                        continue

                    if True in np.isnan(tr1.data):
                        print("Trace contains nan\n",file=output_file)
                        continue

                    if True in np.isnan(tr2.data):
                        print("Trace contains nan\n",file=output_file)
                        continue

                    if True in np.isinf(tr1.data):
                        print("Trace contains inf\n",file=output_file)
                        continue

                    if True in np.isinf(tr2.data):
                        print("Trace contains inf\n",file=output_file)
                        continue


                    # - correlate
                    correlation = cross_covar(tr1.data,tr2.data,
                        max_lag_samples,self.cfg.corr_normalize)[0]
                    
                    # - add to stack
                    if len(correlation) == 2 * max_lag_samples + 1:
                        self._correlations[cp_name]._add_corr(correlation,t)
                    else:
                        print('Empty window or all values zero in window.',
                            file=output_file)
                    
                        




            # - update time
            t += self.cfg.time_window_length - self.cfg.time_overlap

        # - Write results

        for corr in self._correlations.values():
            
            corr.write_stack(output_format=self.cfg.format_output)

        print('Finished a correlation block.')

    def rotate(self,str1,str2,baz1,baz2):

        
        if self.cfg.time_overlap == 0:
            s_temp1 = str1.copy()
            s_temp2 = str2.copy()
        else:
            # deepcopy was taken during slicing already
            s_temp1 = str1
            s_temp2 = str2

        try:
            s_temp1.rotate('NE->RT',back_azimuth=baz1)
        
            
        except:
            print('** data not rotated for stream: ')
            print(s_temp1)
            pass

        try:
            s_temp2.rotate('NE->RT',back_azimuth=baz2)
        except:
            print('** data not rotated for stream: ')
            print(s_temp2)
            pass

        return(s_temp1,s_temp2)



    def update_data(self,t,win_len):


        for trace in self.data:
            
            # is the trace long enough?
            if trace.stats.endtime < (t + win_len):
                
                if self.inv[trace.id] == []:
                    try:
                        self.data.remove(trace)
                        continue
                    except:
                        continue


                # which trace to read next?
                f = self.inv[trace.id].pop(0)
                
                try:
                    
                    # read new trace
                    newtrace = read(f)

                    # add new trace to stream
                    self.data += newtrace
                    self.data._cleanup()
                except:
                    print('** Could not read trace: %s' %trace.id)


                
                # Get the last bit of the trace that we still need
                self.data.trim(starttime=t)

                #self.data._cleanup()
                #gc.collect()
                # Try to create a new stream, to free up memory more successfully:
                #self.data = Stream(traces=self.data.traces)

                

            else:
                pass

            
                # Only trim -- use no more memory than necessary
                #self.data.trim(starttime=t)
        
        # Try to create a new stream, to free up memory more successfully:
        #self.data = Stream(traces=self.data.traces)


    def initialize_data(self):

        self.data = Stream()
        
        

        for channel in self.channels:
            
            
            f = self.inv[channel].pop(0)
            try:
                self.data += read(f)
                
            except IOError:
                print('** problems reading file %s' 
                %self.inv.data[channel])



    def get_prepstring(self):

        prepstring = ''
        if self.cfg.bandpass: 
            prepstring +='b'
        else:
            prepstring += '-'
        if self.cfg.cap_glitch: 
            prepstring += 'g'
        else:
            prepstring += '-'    
        if self.cfg.whiten: 
            prepstring += 'w'
        else:
            prepstring += '-'
        if self.cfg.onebit: 
            prepstring += 'o'
        else:
            prepstring += '-'
        if self.cfg.ram_norm: 
            prepstring += 'r'
        else:
            prepstring += '-'

        return prepstring
        
        

