# Re-vamped correlation script
from __future__ import print_function
from mpi4py import MPI

import os
import time
from ants_2.config import ConfigCorrelation
cfg = ConfigCorrelation()
from ants_2.classes.corrblock import CorrBlock

import ants_2.tools.bookkeep as bk
from obspy import UTCDateTime
from glob import glob
from copy import deepcopy
# 'main':

# - import modules

# - determine own rank, size

# - initialize directories and report text files

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("Hello from rank %g" %rank)
print("Size is %g" %size)





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


    # - get list of files available;
    # - get blocks of channel pairs

    favail = bk.file_avail(cfg)


# - LOOP over blocks:
    
    for b in favail.blocks[rank::size]:

        favail_block = deepcopy(favail)
        # initialize a block of correlations
        c = CorrBlock(b,favail_block,cfg)

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