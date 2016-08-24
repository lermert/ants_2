from __future__ import print_function
import sys
import os
import click
import scripts
from ants_2.config import ConfigDownload, ConfigPreprocess, ConfigCorrelation


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
@click.argument('measure_type',help='Type of measurement.')
@click.option('--speed',help='approx. wave speed')
@click.option('--window',help='window type',default='hann')
@click.option('--plot',default=False)



def run_measurement():

    
    mtype = measr_config['mtype']
    
    # TODo all available misfits --  what parameters do they need (if any.)
    if measr_config['mtype'] in ['ln_energy_ratio','energy_diff']:
        

        g_speed                         =    measr_config['g_speed']
        window_params                   =    {}
        window_params['hw']             =    measr_config['window_params_hw']
        window_params['sep_noise']      =    measr_config['window_params_sep_noise']
        window_params['win_overlap']    =    measr_config['window_params_win_overlap']
        window_params['wtype']          =    measr_config['window_params_wtype']
        window_params['causal_side']    =    measr_config['window_params_causal']
        window_params['plot']           =    measr_config['window_plot_measurements']
    
    measurement(source_config,mtype,step,g_speed=g_speed,window_params=window_params)
    