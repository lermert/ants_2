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
    


@run.command(help='Print directory where ants code is located.')
def show_root():
    from . import _ROOT
    print(_ROOT)


#==============================================================================
# Setup of the directory structure
#==============================================================================

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
    
    
    os.system('cp -r {} ants_code/'.format(_ROOT))
    

    config = ConfigDownload()
    config.initialize()
    config = ConfigPreprocess()
    config.initialize()
    onfig = ConfigCorrelation()
    config.initialize()
    
    
    # .ToDO ... the other configurations follow here...
    
    click.secho('All the input files were copied to ./input, please edit!',color='g')
   
#==============================================================================
# Data download
#==============================================================================

@run.command(help='Download data, input file\
input/config_download.json')
def download():
    from scripts.ant_download import ant_download
    ant_download()

#==============================================================================
# Preprocessing
#==============================================================================

@run.command(help='Remove instrument response, input file input/config_preprocess.json')
def preprocess():
    from scripts.ant_preprocess import preprocess
    preprocess()
    
#==============================================================================
# Correlation
#==============================================================================

@run.command(help='Correlation, input file input/config_correlation.json')
def correlation():
    from scripts.ant_correlation import correlate
    correlate()




#==============================================================================
# Measurement
#==============================================================================

@run.command(help='Take a measurement on the data')
@click.argument('measure_type')
@click.option('--bandpass',help='Filter before measurement',default=None)
@click.option('--speed',help='approx. wave speed in m/s')
@click.option('--hw',help='window half width in seconds')
@click.option('--window',help='window type',default='hann')
@click.option('--plot',help='show plots of the correlations',is_flag=True)
@click.option('--causal',help='In case of using energy_diff measurement, this option selects causal or acausal branch of the correlation',default=True)
@click.option('--sep_noise',help='separation of noise window behind signal window as a multiple of the chosen window halfwidth',default=1)
@click.option('--overlap',help='If set to True, measurements will be taken also if causal and acausal window overlap',default=False)
def measure(measure_type,bandpass,speed,hw,window,plot,causal,sep_noise,overlap):
    from scripts.ant_measurement import measurement
    measure_types = ['ln_energy_ratio','energy_diff']

    if measure_type not in measure_types:
        print('Unrecognized measure_type. measure_type can be: ')
        for t in measure_types:
            print(t)
        return()


    if isinstance(bandpass,unicode):

        bandpass = [float(f) for f in bandpass.split(',')]
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
    
    measurement(measure_type,filt=bandpass,g_speed=g_speed,window_params=window_params)



#==============================================================================
# Determining a source map from a csv file
#==============================================================================


@run.command(help='Plot measured log energy ratios on a map using ray-theoretical kernels')
@click.option('--f',help='Central frequency in Hz')
@click.option('--speed',help='approx. wave speed in m/s')
@click.option('--q',help='Quality factor of surface wave dispersion')
@click.option('--ray_step',help='discretizing step of ray in m',default=1e5)
@click.option('--bin_size',help='Geographic bin size in degree',default=5.)
def sourcemap(f,speed,q,ray_step,bin_size):


    ray_step = float(ray_step) / 1000.
    speed = float(speed)
    f = float(f)
    q = float(q)
    bin_size = float(bin_size)

    if not os.path.exists('data/ln_energy_ratio.measurement.csv'):
        msg = 'Source maps can currently only be created from file ln_energy_ratio.csv. Run ants measure ln_energy_ratio to obtain this file.'
        raise NotImplementedError(msg)

    from ants_2.scripts.ant_sourceimaging import sourcemap
    s = sourcemap('data/ln_energy_ratio.measurement.csv',speed,f,q,ray_step)
    s._temp_kernels()
    s._bin_kernels(bin_size,bin_size)
    s.plot_sourcemap()



#==============================================================================
# Plotting
#==============================================================================

@run.group()
def plot():
    pass
    
@plot.command(help='Plot station map')
@click.option('--bluemarble',help='Plot map background', is_flag=True)
@click.option('--proj',help='Selects matplotlib projection',default='merc',type=str)
def station_map(proj,bluemarble):

    from tools.plot import plot_stations
    plot_stations(projection=proj,bluemarble=bluemarble)



@plot.command(help='Plot stacking of correlation windows')
@click.argument('input_file')
@click.option('--bandpass',help='Bandpass filter: Specify as fmin,fmax,order',default=None,type=str)
@click.option('--pause',help='For animated plots: set pause between frames (in seconds)',default=0.,type=float)
def windows_stack(input_file,bandpass,pause):

    if bandpass is not None:
       try:
           bandpass = [float(nr) for nr in bandpass.split(',')]
       except:
           print('Bandpass format must be: freqmin,freqmax,order')
           bandpass = None
    from tools.plot import plot_converging_stack
    plot_converging_stack(input_file,bandpass,float(pause))



@plot.command(help='Plot single correlation trace')
@click.argument('input_file')
@click.option('--bandpass',help='Bandpass filter: Specify as fmin,fmax,order',default=None,type=str)

def correlation(input_file,bandpass):
    if bandpass is not None:
       try:
           bandpass = [float(nr) for nr in bandpass.split(',')]
       except:
           print('Bandpass format must be: freqmin,freqmax,order')
           bandpass = None
    from tools.plot import plot_correlation
    plot_correlation(input_file,bandpass)




@plot.command(help='Plot section of correlation traces')
@click.argument('directory')
@click.option('--bandpass',help='Bandpass filter: Specify as fmin,fmax,order',default=None,type=str)

def section(directory,bandpass):
    if bandpass is not None:
       try:
           bandpass = [float(nr) for nr in bandpass.split(',')]
       except:
           print('Bandpass format must be: freqmin,freqmax,order')
           bandpass = None
    from tools.plot import plot_section
    plot_section(directory,bandpass)





    