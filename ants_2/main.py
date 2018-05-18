from __future__ import print_function
import sys
import os
import click
import ants_2.scripts
from ants_2.config import ConfigDownload, ConfigPreprocess, ConfigCorrelation
import warnings

warnings.filterwarnings("ignore",message='Found more than one matching coordinates. Returning first.')

@click.group()
def run():
    """
    Main routine for noise correlation modeling and noise source inversion.
    """
    pass
    


# @run.command(help='Print directory where ants code is located.')
# def show_root():
#     from . import _ROOT
#     print(_ROOT)


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
    os.mkdir(os.path.join('.','data','treated'))
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
    config = ConfigCorrelation()
    config.initialize(check_params=False)
    
    
    # .ToDO ... the other configurations follow here...
    
    click.secho('All the input files were copied to ./input, please edit!',color='g')
   
#==============================================================================
# Data download
#==============================================================================

@run.command(help='Download data\ninput file\
input/config_download.json')
def download():
    from ants_2.scripts.ant_download import ant_download
    ant_download()

#==============================================================================
# Preprocessing
#==============================================================================

@run.command(help='Remove instrument response\ninput file input/config_preprocess.json')
def preprocess():
    from ants_2.scripts.ant_preprocess import preprocess
    preprocess()
    
#==============================================================================
# Correlation
#==============================================================================

@run.command(help='Correlation\ninput file input/config_correlation.json')
def correlation():
    from ants_2.scripts.ant_correlation import correlate
    correlate()




#==============================================================================
# Measurements
#==============================================================================

@run.group(help='Take a measurement on the data.')
def measure():
    pass


@measure.command()
@click.option('--bandpass',help='freqmin,freqmax,order: Butterworth bandpass filter',default=None)
@click.option('--speed',help='approx. wave speed in m/s',default=None)
@click.option('--hw',help='window half width in seconds',default=None)
@click.option('--window',help='window type',default='hann')
@click.option('--plot',help='show plots of the correlations',is_flag=True)
@click.option('--sep_noise',help='separation of noise window behind signal window as a multiple of the chosen window halfwidth',default=1)
@click.option('--overlap',help='If set to True, measurements will be taken even if causal and acausal window overlap',default=False)

def ln_energy_ratio(bandpass,speed,hw,window,plot,sep_noise,overlap):
    """

    Measure logarithmic energy ratio.
    
    Example: ants measure --speed 3000.0 --hw 100 --bandpass 0.1,0.2,4
    """

    from ants_2.scripts.ant_measurement import measurement
    

    if isinstance(bandpass,unicode):
        try:
            bandpass = [float(f) for f in bandpass.split(',')]
            print("Filtering between %g and %g Hz." %(bandpass[0],bandpass[1]))
        except:
           print('Bandpass format must be: freqmin,freqmax,order')
           bandpass = None
    
    
    # TODo all available misfits --  what parameters do they need (if any.)
    
        
    if speed is not None:
        speed                           =    float(speed)

    window_params                       =    {}

    if hw is not None:
        window_params['hw']             =    float(hw)
    else:
        window_params['hw']             =    hw

    window_params['sep_noise']          =    int(sep_noise)
    window_params['win_overlap']        =    bool(overlap)
    window_params['wtype']              =    window
    window_params['plot']               =    bool(plot)
    window_params['causal_side']        =    True
    
    measurement(mtype='ln_energy_ratio',filt=bandpass,
        g_speed=speed,window_params=window_params)

#@click.option('--causal',help='In case of using energy_diff measurement, this option selects causal or acausal branch of the correlation',default=True)

#==============================================================================
# Determining a source map from a csv file
#==============================================================================


@run.command()
@click.option('--f',help='Central frequency in Hz')
@click.option('--speed',help='approx. wave speed in m/s')
@click.option('--q',help='Quality factor of surface wave dispersion')

@click.option('--ray_step',help='discretizing step of ray in m',default=1e5)
@click.option('--bin_size',help='Geographic bin size in degree',default=5.)
@click.option('--min_snr',help='Minimum signal to noise ratio for meausurements',default=0.0)
@click.option('--box',help='Geographic limits for map: latmin,latmax,lonmin,lonmax in degree',default=None)
def sourcemap(f,speed,q,ray_step,bin_size,min_snr,box):
    
    """
    Plot measured log energy ratios on a map using ray-theoretical kernels.

    Example: ants sourcemap --f 0.14 --speed 3000 --q 120
    """


    if not os.path.exists('data/ln_energy_ratio.measurement.csv'):
        msg = 'Source maps can currently only be created from file ln_energy_ratio.csv.\nRun ants measure ln_energy_ratio to obtain this file.'
        raise NotImplementedError(msg)


    try:
        
        speed = float(speed)
        f = float(f)
        q = float(q)
    except:
        msg = 'For this measurement, options f, speed and q must be specified.'
        raise ValueError(msg)
    
    bin_size = float(bin_size)
    ray_step = float(ray_step) / 1000.
    min_snr = float(min_snr)

    if isinstance(box,unicode):
        try:
            box = [float(l) for l in box.split(',')]
            print("Geographic limits: {} {} {} {}".format(*box))
        except:
           print('Geographic limits format: latmin,latmax,lonmin,lonmax.\
        lat must be in (-90,90) and lon in (-180,180)')
           box = [-90.,89.,-180,179]
    else:
        box = [-80.,80.,-179,179]

   
    latmin = box[0]
    latmax = box[1]
    lonmin = box[2]
    lonmax = box[3]

    if lonmin > lonmax:
        raise ValueError('Geographic limits format: latmin,latmax,lonmin,lonmax.\
        lat must be in (-90,90) and lon in (-180,180)')

    from ants_2.scripts.ant_sourceimaging import sourcemap
    s = sourcemap('data/ln_energy_ratio.measurement.csv',speed,f,q,ray_step,min_snr)
    s._temp_kernels()
    s._bin_kernels(bin_size,bin_size,lonmin,lonmax,latmin,latmax)
    s.plot_sourcemap()



#==============================================================================
# Plotting
#==============================================================================

@run.group(help='Various plotting functions')
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
@click.option('--bandpass',help='freqmin,freqmax,order: Butterworth bandpass filter',default=None,type=str)
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
@click.option('--bandpass',help='freqmin,freqmax,order: Butterworth bandpass filter',default=None,type=str)

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
@click.option('--bandpass',help='freqmin,freqmax,order: Butterworth bandpass filter',default=None,type=str)

def section(directory,bandpass):
    if bandpass is not None:
       try:
           bandpass = [float(nr) for nr in bandpass.split(',')]
       except:
           print('Bandpass format must be: freqmin,freqmax,order')
           bandpass = None
    from tools.plot import plot_section
    plot_section(directory,bandpass)





    
