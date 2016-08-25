import os
import numpy as np
import pandas as pd
from math import log, pi
import click
import json
from scipy.signal import hilbert
from glob import glob
from obspy import read, Trace
from obspy.geodetics import gps2dist_azimuth
import matplotlib.pyplot as plt
from ants_2.tools import measurements as rm
from ants_2.tools.windows import get_window, my_centered, snratio



def get_station_info(stats):

    sta1 = '{}.{}.{}.{}'.format(stats.network,stats.station,stats.location,
    stats.channel)
    sta2 = '{}.{}.{}.{}'.format(stats.sac.kuser0.strip(),stats.sac.kevnm.strip(),
    stats.sac.kuser1.strip(),stats.sac.kuser2.strip())
    lat1 = stats.sac.stla
    lon1 = stats.sac.stlo
    lat2 = stats.sac.evla
    lon2 = stats.sac.evlo
    dist = stats.sac.dist
    az = gps2dist_azimuth(lat1,lon1,lat2,lon2)[2]
    
    
    return([sta1,sta2,lat1,lon1,lat2,lon2,dist,az])

     

def measurement(mtype,**options):
    
    """
    Get measurements on noise correlation data and synthetics. 
    options: g_speed,window_params (only needed if mtype is ln_energy_ratio or enery_diff)
    """
    
    
    files = glob(os.path.join('data','correlations','*.SAC'))
    
    
    columns = ['sta1','sta2','lat1','lon1','lat2','lon2','dist','az',
    'obs','snr']
    measurements = pd.DataFrame(columns=columns)
    
    
    
    if files == []:
        msg = 'No input found!'
        raise ValueError(msg)
    
    i = 0
    with click.progressbar(files,label='Taking measurements...') as bar:
        
        for f in bar:
            
            try: 
                tr_o = read(f)[0]

            except:
                print('\nCould not read data: '+os.path.basename(f))
                i+=1
                continue
           
            
            # Get all the necessary information
            info = get_station_info(tr_o.stats)
           
            # Take the measurement
           
            func = rm.get_measure_func(mtype)
            msr_o = func(tr_o,**options)
            msr = msr_o
            snr = snratio(tr_o,**options)
            
            
            info.extend([msr,snr])
            measurements.loc[i] = info

            # step index
            i+=1
    
    filename = '{}.measurement.csv'.format(mtype)
    measurements.to_csv(os.path.join('data',filename),index=None)

