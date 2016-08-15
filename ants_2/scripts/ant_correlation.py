# Re-vamped correlation script
from __future__ import print_function
from mpi4py import MPI

import os

from ants_2.config import ConfigCorrelation
cfg = ConfigCorrelation()

from obspy import UTCDateTime
from glob import glob
# 'main':

# - import modules

# - determine own rank, size

# - initialize directories and report text files

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("Hello from rank %g" %rank)
print("Size is %g" %size)



class file_inventory(object):

    """
    For each station id, keep a dictionary of available files within the time range requested in cfg.
    """

    def __init__(self,cfg):

        
        indirs = cfg.indirs
        t0 = UTCDateTime(cfg.time_start) 
        t1 = UTCDateTime(cfg.time_end)
        fileformat = cfg.input_format

        self.get_processed_data(indirs,t0,t1,fileformat)


    def get_processed_data(self,indirs,t0,t1,filefmt):
        
        files = []
        self.data = {}
        

        for indir in indirs:
            files.extend(glob(os.path.join(indir,'*.'+filefmt.lower())))
            files.extend(glob(os.path.join(indir,'*.'+filefmt.upper())))
            
        
        for f in files:
            print(f)
            fn = os.path.basename(f).split('.')
            st = UTCDateTime('{}-{}T{}:{}:{}'.format(*fn[4:9]))
            et = UTCDateTime('{}-{}T{}:{}:{}'.format(*fn[9:14]))
            print(et)
            if st > t1 or et < t0:
                continue
            else:

                station = '{}.{}'.format(*fn[0:2])
                channel = '{}.{}.{}.{}'.format(*fn[0:4])

                if station not in self.data.keys():
                    self.data[station] = {}

                if channel not in self.data[station].keys():
                    self.data[station][channel] = []

                self.data[station][channel].append(f)

        return()

    def get_blocks(self,n):

        # ToDo check if inventory already there?
        staids = self.data.keys()
        staids.sort()
        # ToDo build in updating mode again?
        blcks = []
        idprs = []

        n_ids = len(staids)
        
        for i in range(n_ids):
            for j in range(i+1,n_ids):

                if len(idprs) == n:
                    blcks.append(idprs)
                    idprs = []

                idprs.append((staids[i],staids[j]))

        if len(idprs) < n:
            blcks.append(idprs)

        return blcks


def correlate():
  # Create output directory, if necessary

    outdir = os.path.join('data','processed')
     
    if rank == 0 and not os.path.exists(outdir):
        os.mkdir(outdir)

    comm.Barrier()
    # Create own output directory, if necessary

    rankdir = os.path.join(outdir,'rank_%g' %rank)
    if not os.path.exists(rankdir):
        os.mkdir(rankdir)


    # correlation report file

    output_file = os.path.join(rankdir,
        'correlation_report_rank%g.txt' %rank)

    if os.path.exists(output_file):
        ofid = open(output_file,'a')
        print('UPDATING, Date:',file=ofid)
        print(time.strftime('%Y.%m.%dT%H:%M'),file=ofid)
    else:
        ofid = open(output_file,'w')
        print('CORRELATING, Date:',file=ofid)
        print(time.strftime('%Y.%m.%dT%H:%M'),file=ofid)


    # - get list of files available

    content = file_inventory(cfg)

    # - determine station pairs;

    blocks  = content.get_blocks()

# - LOOP over blocks:

# - block.correlate
 
# - block.write

# - move the calculated correlations to the output directory

#class fileentry(object):

 #   def __init__(self,filepath):
 #       self.filepath()



    # List all the files 

    # loop through files:
    # - if file ends before t0 or starts after t1, continue
    # - register channel; all files of one channel belong to this channel (dictionary?)

    # form a set of channels
    # group channels by station
    #
    # How to attach time range of that file?
    # Could be useful for plotting availability.