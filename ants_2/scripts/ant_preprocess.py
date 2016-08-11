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
        print('UPDATING, Date:',file=ofid)
        print(time.strftime('%Y.%m.%dT%H:%M'),file=ofid)
    else:
        ofid = open(output_file,'w')
        print('PROCESSING, Date:',file=ofid)
        print(time.strftime('%Y.%m.%dT%H:%M'),file=ofid)


    # select input files for this rank    
    content = content[rank::size]
    if cfg.testrun: # Only 3 files randomly selected
        indices = randint(0,len(content),3)
        content = [content[j] for j in indices]

    # Loop over input files
    for filepath in content:

        print('-------------------------------------',file=ofid)
        print('Attempting to process:',file=ofid)
        print(os.path.basename(filepath),file=ofid)
        
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
            
    