# A script to process ambient vibration records
from __future__ import print_function
# Use the print function to be able to switch easily between stdout and a file
from mpi4py import MPI
import os
import sys
import time

#from obspy import read, Stream,  Trace, UTCDateTime
#from glob import glob
from numpy.random import randint
from ants_2.tools.bookkeep import find_raw_files
from ants_2.config import ConfigPreprocess
cfg = ConfigPreprocess()
from ants_2.classes.prepstream import PrepStream


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("Hello from rank %g" %rank)
print("Size is %g" %size)


def preprocess():
    """
    
    This script preprocesses the MSEED files in the input directories 
    specified in the input file.
 
    
    """


    # Create output directory, if necessary

    outdir = os.path.join('data','processed')
     
    if rank == 0 and not os.path.exists(outdir):
        os.mkdir(outdir)
    
    comm.Barrier()
    
       

    # Create own output directory, if necessary

    rankdir = os.path.join(outdir,
        'rank_%g' %rank)
    if not os.path.exists(rankdir):
        os.mkdir(rankdir)

    
    #- Find input files
    
    content = find_raw_files(cfg.input_dirs,
        cfg.input_format)
       


    # processing report file

    output_file = os.path.join(rankdir,
        'processing_report_rank%g.txt' %rank)

    if os.path.exists(output_file):
        ofid = open(output_file,'a')
        print('PROCESSING, Date:',file=ofid)
        print(time.strftime('%Y.%m.%dT%H:%M'),file=ofid)
    else:
        ofid = open(output_file,'w')
        print('UPDATING, Date:',file=ofid)
        print(time.strftime('%Y.%m.%dT%H:%M'),file=ofid)


    # select input files for this rank    
    content = content[rank::size]
    if cfg.testrun: # Only 3 files randomly selected
        indices = randint(0,len(content),3)
        content = [content[j] for j in indices]

    # Loop over input files
    for filepath in content:
        
        try:
            prstr = PrepStream(filepath,ofid)
        except:
            print('** Problem opening file, skipping: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue

        if len(prstr.stream) == 0:
            print('** No data in file, skipping: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue
        
        try:
            prstr.prepare(cfg)
        except:
            print('** Problems preparing stream: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue
            
        try:
            prstr.process(cfg)
        except:
            print('** Problems processing stream: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue

        try:
            prstr.write(rankdir,cfg)
        except:
            print('** Problems writing stream: ',file=ofid)
            print('** %s' %filepath,file=ofid)

        
    ofid.close()

    print("Rank %g has completed processing." 
        %rank,file=None)
    
    
    try:
        os.system('mv '+rankdir+'/* '+outdir)
    except:
        pass

    os.system('rmdir '+rankdir)
            
    # def getfilepath(self,rankdir,stats,startonly=False):
        
    #     inf = [
    #         stats.network,
    #     stats.station,
    #     stats.location,
    #     stats.channel,
    #     stats._format
    #     ]
        
            
    #     t1=stats.starttime.strftime('%Y.%j.%H.%M.%S')
    #     t2=stats.endtime.strftime('%Y.%j.%H.%M.%S')
    #     if startonly: 
    #         t2 = '*'
        
    #     inf.append(t1,t2)
        
    #     #yr=str(t1[0:4])
        
    #     filenew = '{}.{}.{}.{}.{}.{}.{}'.format(inf)
    #     filepathnew = os.path.join(rankdir,filenew)
    #     #filepathnew=rankdir+'/'+network+'.'+station+'.'+location+'.'+\
    #     #channel+'.' + t1 + '.' +t2+'.'+prepname+'.'+format 
        
    #     return filepathnew
            
    
    # def find_rawdatafiles(self):
        
    #     indirs = cfg.input_dirs
    #     format = cfg.input_format.lower()
    #     content=list()
    #     for indir in indirs:
    #         print(indir)
    #         content.extend(glob(os.path.join(indir,'*'+format)))
    #         content.extend(glob(os.path.join(indir,'*'+format.upper())))
            
    #     content.sort()
    #     print(content)
    #     return content
        
    # def read_file(self,filepath,ofid):
        
    #     if cfg.verbose:
    #         print('===========================================================',\
    #         file=ofid)
    #         print('* opening file: '+filepath+'\n',file=ofid)
            
    #     #- read data
    #     try:
    #         data=read(filepath)
    #     except (TypeError, IOError):
    #         if cfg.verbose==True: print('** file wrong type or not found, skip.',
    #         file=ofid)
    #         return []
    #     except:
    #         if cfg.verbose: print('** unexpected read error, skip.',
    #         file=ofid)
    #         return []
            
    #     #- check if this file contains data
    #     if len(data) == 0:
    #         print('** file contains no data, skip.',file=ofid)
    #         return []
        
    #     return data
    
    
    
    # 