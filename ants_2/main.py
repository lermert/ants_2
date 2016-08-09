from __future__ import print_function
import sys
import os
import click
import scripts
from ants_2.config import ConfigDownload, ConfigPreprocess


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
    
    
    # .ToDO ... the other configurations follow here...
    
    click.secho('All the input files were copied to ./input, please edit!',color='g')
   
@run.command(help='Download data from IRIS or Arclink or both.\nEdit input file\
input/config_download.json')
def par_download():
    from scripts.ant_download import ant_download
    ant_download()

@run.command(help='Instrument correction and other preprocessing.')
def par_preprocess():
    from scripts.ant_preprocess import preprocess
    preprocess()
    