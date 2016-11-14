import os
import numpy as np
import pandas as pd
from math import log, pi, isnan
import click
import json
from scipy.signal import hilbert
from glob import glob
from obspy import read, Trace
from obspy.geodetics import gps2dist_azimuth
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
    az,baz = gps2dist_azimuth(lat1,lon1,lat2,lon2)[1:]
    
    
    return([sta1,sta2,lat1,lon1,lat2,lon2,dist,az,baz])

     
# ToDo: Channel choice
def measurement(mtype,filt,dir,cha1='',cha2='',**options):
    
    """
    Get measurements on noise correlation data and synthetics. 
    options: g_speed,window_params (only needed if mtype is ln_energy_ratio or enery_diff)
    """
    
    
    files = glob(os.path.join(dir,'*{}*{}*.SAC'.format(cha1,cha2)))
    
    
    columns = ['sta1','sta2','lat1','lon1','lat2','lon2','dist','az',
    'baz','obs','snr']
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
                continue
           
            # Filter
            if filt is not None:
                tr_o.taper(type='cosine',max_percentage=0.05)
                tr_o.filter('bandpass',freqmin=filt[0],freqmax=filt[1],
                    corners=filt[2],zerophase=True)
            
            # Get all the necessary information
            info = get_station_info(tr_o.stats)

           
            # Take the measurement
           
            func = rm.get_measure_func(mtype)
            msr_o = func(tr_o,**options)
            msr = msr_o
            snr = snratio(tr_o,**options)

            if isnan(msr):
                continue
            else:
                info.extend([msr,snr])
                measurements.loc[i] = info

                # step index
                i+=1
    
    filename = '{}.measurement.csv'.format(mtype)
    measurements.to_csv(os.path.join('data',filename),index=None)

    fh = open(os.path.join('data','measurement_options.txt'),'w')

    if filt is not None:
        fh.write("Butterworth filter freqmin,freqmax,order: ")
        fh.write("{},{},{}\n".format(*filt))
    for key,value in options.items():
        fh.write("{}: {}\n".format(key,value))

    fh.close()



