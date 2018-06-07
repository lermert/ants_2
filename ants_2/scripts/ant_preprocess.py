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
from obspy import UTCDateTime
from obspy.clients.fdsn.client import Client
from obspy.core.event.catalog import Catalog
from ants_2.tools.bookkeep import find_files
from ants_2.tools.prepare import get_event_filter
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
    if rank == 0 and cfg.verbose:
        print(cfg.__dict__)
    
    comm.Barrier()

    event_filter = None

    if cfg.gcmt_exclude:

        if rank == 0:
            c = Client()
            cata = c.get_events(starttime=UTCDateTime(cfg.gcmt_begin),
                endtime=UTCDateTime(cfg.gcmt_end),catalog='GCMT',
                minmagnitude=5.6)
    
            event_filter = get_event_filter(cata,cfg.Fs_new[-1],
                t0=UTCDateTime(cfg.gcmt_begin),
                t1=UTCDateTime(cfg.gcmt_end))

        
        # communicate event_filter (would it be better 
        # if every rank sets it up individually?)
        event_filter = comm.bcast(event_filter,root=0)
    
    if cfg.event_exclude_local_cat:

        local_cat = Catalog()
        
        if rank == 0:
            c = Client()
            local_cat.extend(c.get_events(
                    starttime=UTCDateTime(cfg.event_exclude_local_cat_begin),
                    endtime=UTCDateTime(cfg.event_exclude_local_cat_end),
                    #catalog=catalog,
                    minmagnitude=cfg.event_exclude_local_cat_minmag,
                    latitude=cfg.event_exclude_local_cat_lat,
                    longitude=cfg.event_exclude_local_cat_lon,
                    maxradius=cfg.event_exclude_local_cat_radius))
            print(len(local_cat),"events in local earthquake catalog.")
        # communicate event_filter (would it be better 
        # if every rank sets it up individually?)
        local_cat = comm.bcast(local_cat,root=0)

    # Create own output directory, if necessary
    rankdir = os.path.join(outdir,
        'rank_%g' %rank)
    if not os.path.exists(rankdir):
        os.mkdir(rankdir)

    
    #- Find input files
    
    content = find_files(cfg.input_dirs,
        cfg.input_format)
    if rank==0:
        print(len(content), "files found") 
    #print(content)

    # processing report file
    sys.stdout.flush()
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
            prstr.process(cfg,event_filter,local_cat)
        except:
            print('** Problems processing stream: ',file=ofid)
            print('** %s' %filepath,file=ofid)
            continue

        try:
            prstr.write(rankdir,cfg)
        except:
            print('** Problems writing stream: ',file=ofid)
            print('** %s' %filepath,file=ofid)

        ofid.flush()
        
    ofid.close()

    print("Rank %g has completed processing." 
        %rank,file=None)
    
    
    try:
        os.system('mv '+rankdir+'/* '+outdir)
    except:
        pass

    os.system('rmdir '+rankdir)
            
    
