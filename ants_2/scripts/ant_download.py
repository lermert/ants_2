#!/usr/bin/env python

# A script to download ambient vibration records
from __future__ import print_function
import os
import sys
import shutil
from obspy import UTCDateTime
from math import ceil
import json

from mpi4py import MPI
from glob import glob
from ants_2.config import ConfigDownload

cfg = ConfigDownload()
from .. import _ROOT
from obspy.clients import fdsn


def ant_download():
 #===============================================================================
        # preliminaries
        #===============================================================================
        
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    outdir = os.path.join('data','raw')
    targetloc=os.path.join(outdir,'rank'+str(rank))
    if not os.path.isdir(targetloc):
        os.mkdir(targetloc)
    respfileloc=os.path.join('meta','resp')
    if os.path.isdir(respfileloc)==False:
        cmd='mkdir '+respfileloc
        os.system(cmd)
    
    if rank == 0:
        client = fdsn.Client()
         #===============================================================================
        #- read station list
        #- create output directory
        #- set parameters #===============================================================================
 
    # network, channel, location and station list
    #stalist=cfg.ids#os.path.join('input','downloadlist.txt')
    fh=open(cfg.ids,'r')
    ids=fh.read().split('\n')
    
    # Verbose?
    if cfg.verbose:
        v=True
        vfetchdata='-v '
    else:
        vfetchdata=''
        
    # Quality?
    quality = cfg.quality
        
    # time interval of request
    t1=cfg.t_start
    t1str=UTCDateTime(t1).strftime('%Y.%j.%H.%M.%S')
    t2=cfg.t_end
    t2str=UTCDateTime(t2).strftime('%Y.%j.%H.%M.%S')
    
    # data segment length
    if cfg.seconds_segment==None:
        winlen=UTCDateTime(t2)-UTCDateTime(t1)
    else:
        winlen = int(cfg.seconds_segment)
    
    # minimum length
    minlen = int(cfg.seconds_minimum)
    
    # geographical region
    lat_min=cfg.lat_min
    lat_max=cfg.lat_max
    lon_min=cfg.lon_min
    lon_max=cfg.lon_max
    
    #===============================================================================
    
    #- Assign each rank its own chunk to download
    #===============================================================================
    
    clen=int(float(len(ids))/float(size))
    chunk=(rank*clen, (rank+1)*clen)
    myids=ids[chunk[0]:chunk[1]]
    if rank==size-1:
        myids=ids[chunk[0]:]
       
    
    #===============================================================================
    
    # Station loop
    #===============================================================================
      
    for id in myids:
        
        if id=='': continue
        network=id.split('.')[0]
        station=id.split('.')[1]
        channel=id.split('.')[3]
         #===============================================================================
        
        # Time window loop
         #===============================================================================
        t = UTCDateTime(t1)
        if not cfg.download_response_only:
            while t < UTCDateTime(t2):
                
                tstart = UTCDateTime(t).strftime('%Y-%m-%d')
                tstartstr = UTCDateTime(t).strftime('%Y.%j.%H.%M.%S')
                
                tstep = min((UTCDateTime(t)+winlen),UTCDateTime(t2)).\
                strftime('%Y-%m-%d')
                tstepstr = min((UTCDateTime(t)+winlen),UTCDateTime(t2)).\
                strftime('%Y.%j.%H.%M.%S')
                
                
                #-Formulate a polite request
                filename=os.path.join(targetloc,id+'.'+tstartstr+'.'+tstepstr+'.mseed')
                
                  
                if os.path.exists(filename)==False:
                    #print network, station, location, channel
                    print('\n Rank '+str(rank),file=None)
                    print('\n Attempting to download data from: '+id,file=None)
                    print(filename)
                    
                    reqstring_iris = '{} {} -N {} -S {} -C {} -s {} -e {} -msl {} --lat \
                    {}:{} --lon {}:{} -o {} -Q {}'.format(os.path.join(_ROOT,'tools_ext','FetchData')\
                    ,vfetchdata,network,station,channel,tstart,tstep,minlen,lat_min,lat_max,lon_min,\
                    lon_max,filename,quality)
                    
                    reqstring_arclink = '{} {} -N {} -S {} -C {} -s {} -e {} -msl {} --lat \
                    {}:{} --lon {}:{} -o {} -Q {}'.format(os.path.join(_ROOT,'tools_ext','FetchDataArc')\
                    ,vfetchdata,network,station,channel,tstart,tstep,minlen,lat_min,lat_max,lon_min,\
                    lon_max,filename,quality)
                    
                    
                    #reqstring=_ROOT+'/tools/FetchData '+vfetchdata+' -N '+network+ \
                    # ' -S '+station+' -C '+channel+' -s '+tstart+' -e '+tstep+ \
                    # ' -msl '+minlen+' --lat '+lat_min+':'+lat_max+ \
                    #' --lon '+lon_min+':'+lon_max+' -o '+filename+' -Q '+quality
                    if cfg.data_center == 'iris' or cfg.data_center=='any':
                        os.system(reqstring_iris)
                    elif cfg.data_center == 'arclink' or cfg.data_center=='any': 
                        os.system(reqstring_arclink)
                t += winlen
        
        tstart = UTCDateTime(t1).strftime('%Y-%m-%d')
        tstep = UTCDateTime(t2).strftime('%Y-%m-%d')
        print('\n Downloading response information from: '+id+'\n')
        
        #===============================================================================
        
        # Within Station loop: Download resp files
         #===============================================================================        
        reqstring_resp_iris = '{} {} -N {} -S {} -C {} -s {} -e {} --lat \
        {}:{} --lon {}:{} -rd {} -Q {}'.format(
        os.path.join(_ROOT,'tools_ext','FetchData'),vfetchdata,network\
        ,station,channel,tstart,tstep,lat_min,lat_max,lon_min,\
        lon_max,respfileloc,quality)
        
        reqstring_resp_arclink = '{} {} -N {} -S {} -C {} -s {} -e {} --lat \
        {}:{} --lon {}:{} -rd {} -Q {}'.format(
        os.path.join(_ROOT,'tools_ext','FetchDataArc'),vfetchdata,network\
        ,station,channel,tstart,tstep,lat_min,lat_max,lon_min,\
        lon_max,respfileloc,quality)
        
        if cfg.data_center == 'iris' or cfg.data_center=='any':
            os.system(reqstring_resp_iris)
        elif cfg.data_center == 'arclink' or cfg.data_center=='any':
            os.system(reqstring_resp_arclink)
            
        

    # Clean up (some files come back with 0 data)
    os.system(os.path.join(_ROOT,'tools','cleandir.sh')+' '+targetloc)
    cmd = 'mv '+targetloc+'/* '+targetloc+'/..'
    print(cmd)
    os.system(cmd)
    os.system('rmdir '+targetloc)

        #===============================================================================
        
        # Separate Station loop: Download stationxml
         #===============================================================================        
    if rank == 0:
        for id in ids:
            if id=='': continue
            network=id.split('.')[0]
            station=id.split('.')[1]
            xmlfile=os.path.join('meta','stationxml','{}.{}.xml'.format(network,station))
            # Metadata request with obspy
            if os.path.exists(xmlfile)==False:
                client.get_stations(network=network,station=station,
                filename=xmlfile,level='response')        

    comm.Barrier()
 #===============================================================================
 
     # After download completed on all ranks: Check availability
  #==============================================================================   
 
    
    if rank==0:
        outfile=os.path.join(outdir,'download_report.txt')
        outf=open(outfile,'w')
        
        print('Attempted to download data from stations: \n',file=outf)
        print('****************************************** \n',file=outf)
        for id in ids:
            print(id,file=outf)
        print('****************************************** \n',file=outf)
        stalist=os.path.join('input','downloadlist.txt')
        fh=open(stalist,'r')
        ids=fh.read().split('\n')
        
        noreturn=[]
        
        for id in ids:
            if id=='': continue
           
            fls=glob(os.path.join(outdir,id+'*'))
            fls.sort()
           
            if fls != []:
                print('Files downloaded for id: '+id,file=outf)
                print('First file: '+fls[0],file=outf)
                print('Last file: '+fls[-1],file=outf)
                print('****************************************** \n',file=outf)    
            else: 
                noreturn.append(id)
            
        if noreturn != []:
            print('NO files downloaded for: \n',file=outf)
            
            print(noreturn,file=outf)
     
        print('****************************************** \n',file=outf)
        print('Download parameters were: \n',file=outf)
        print('****************************************** \n',file=outf)
        outf.close()
        
        os.system('cat input/config_download.json >> '+outfile)
    
    return()
            