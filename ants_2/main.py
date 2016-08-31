from __future__ import print_function
import sys
import os
import click
import scripts
from ants_2.config import ConfigDownload, ConfigPreprocess, ConfigCorrelation
import warnings

warnings.simplefilter("ignore")

@click.group()
def run():
    """
    Main routine for noise correlation modeling and noise source inversion.
    """
    pass
    

@run.command(help='Initialize folder for a new project.')
def new_project():
    if not os.listdir('.')==[]:
        click.echo('Project exists already, can only start in an empty directory.')
        exit()
    from . import _ROOT
    #ToDo put input files at end so we can use the existing if the json doesn't exist
    os.mkdir(os.path.join('.','input'))
    os.mkdir(os.path.join('.','data'))
    os.mkdir(os.path.join('.','data','raw'))
    os.mkdir(os.path.join('.','data','processed'))
    os.mkdir(os.path.join('.','data','correlations'))
    os.mkdir(os.path.join('.','meta'))
    os.mkdir(os.path.join('.','meta','resp'))
    os.mkdir(os.path.join('.','meta','stationxml'))
    os.mkdir(os.path.join('.','ants_code'))
    
    
    os.system('cp -r {} ants_code/'.format(os.path.join(_ROOT,'scripts')))
    os.system('cp -r {} ants_code/'.format(os.path.join(_ROOT,'tools')))
    config = ConfigDownload()
    config.initialize()
    config = ConfigPreprocess()
    config.initialize()
    onfig = ConfigCorrelation()
    config.initialize()
    
    
    # .ToDO ... the other configurations follow here...
    
    click.secho('All the input files were copied to ./input, please edit!',color='g')
   
@run.command(help='Download data from IRIS or Arclink or both.\nEdit input file\
input/config_download.json')
def download():
    from scripts.ant_download import ant_download
    ant_download()

@run.command(help='Instrument correction and other preprocessing.')
def preprocess():
    from scripts.ant_preprocess import preprocess
    preprocess()
    

@run.command(help='Correlation')
def correlation():
    from scripts.ant_correlation import correlate
    correlate()


@run.command(help='Plot stations for which data is available.')
@click.option('--bluemarble', default=False)
@click.option('--proj',default='merc')
def plot_stations(proj,bluemarble):
    from tools.plot import plot_stations
    plot_stations(projection=proj,bluemarble=bluemarble)


@run.command(help='Take a measurement on the data.')
@click.argument('measure_type')
@click.option('--bandpass',help='Filter before measurement',default=None)
@click.option('--speed',help='approx. wave speed in m/s')
@click.option('--hw',help='window half width in seconds')
@click.option('--window',help='window type',default='hann')
@click.option('--plot',default=False)
@click.option('--causal',default=True)
@click.option('--sep_noise',default=1)
@click.option('--overlap',default=False)
def measure(measure_type,bandpass,speed,hw,window,plot,causal,sep_noise,overlap):
    from scripts.ant_measurement import measurement
    measure_types = ['ln_energy_ratio','energy_diff']

    if measure_type not in measure_types:
        print('Unrecognized measure_type. measure_type can be: ')
        for t in measure_types:
            print(t)
        return()


    if isinstance(bandpass,unicode):

        filt = [float(f) for f in bandpass.split(',')]
        print("Filtering between %g and %g Hz." %(filt[0],filt[1]))

    
    
    # TODo all available misfits --  what parameters do they need (if any.)
    if measure_type in ['ln_energy_ratio','energy_diff']:
        # ToDo check whether speed and hw options are passed in.
        

        g_speed                         =    float(speed)
        window_params                   =    {}
        window_params['hw']             =    float(hw)
        window_params['sep_noise']      =    int(sep_noise)
        window_params['win_overlap']    =    bool(overlap)
        window_params['wtype']          =    window
        window_params['causal_side']    =    bool(causal)
        window_params['plot']           =    bool(plot)
    
    measurement(measure_type,filt=filt,g_speed=g_speed,window_params=window_params)
    