import io
import json
import os
import click


DEFAULT_Download = {
    "verbose":True,
    "data_center":"any",
    "data_user":"hugo.bruno.kunz@gmail.com",
    "ids":"input/downloadlist.txt",
    "quality":"?",
    "t_start":"1970-01-01,00:00:00",
    "t_end":"1971-01-01,00:00:00",
    "seconds_segment":None,
    "seconds_minimum":0,
    "lat_min":-90.,
    "lat_max":90.,
    "lon_min":-180.,
    "lon_max":180.
}

CONFIG_Download = os.path.join('input','config_download.json')

class ConfigDownload(object):
    """Contains basic parameters for the job (paths, etc.)"""

    def __init__(self):
        self.verbose = None
        self.data_center = None
        self.data_user = None
        self.ids = None
        self.quality = None
        self.t_start = None
        self.t_end = None
        self.seconds_segment = None
        self.seconds_minimum = None
        self.lat_min = None
        self.lat_max = None
        self.lon_min = None
        self.lon_max = None
        
        self.initialize()
        
    def initialize(self):
        """Populates the class from ./config.json.
        If ./config.json does not exist, writes a default file and exits.
        """

        if not os.path.exists(CONFIG_Download):
            with io.open(CONFIG_Download, 'wb') as fh:
                json.dump(DEFAULT_Download, fh, sort_keys=True, indent=4, separators=(",", ": "))
            return()

        # Load all options.
        with io.open(CONFIG_Download, 'r') as fh:
            data = json.load(fh)
            
        for key, value in data.iteritems():
            setattr(self, key, value)




DEFAULT_Preprocess = {
    "verbose":True,
    "testrun":False,
    "input_dirs":[],
    "input_format":"MSEED",
    "quality_minlengthsec":0.,
    "quality_maxgapsec":0.,
    "event_exclude":False,
    "event_exclude_winsec":[],
    "event_exclude_std":2.,
    "event_exclude_n":4,
    "event_exclude_freqmin":0.01,
    "event_exclude_freqmax":1.0,
    "event_exclude_level":2.,
    "wins":False,
    "wins_len_sec":8192,
    "wins_trim":True,
    "wins_detrend":True,
    "wins_demean":True,
    "wins_taper":0.05,
    "wins_cap":False,
    "wins_cap_threshold":15,
    "Fs_old":[],
    "Fs_new":[],
    "Fs_antialias_factor":0.4,
    "instr_correction":True,
    "instr_correction_unit":'VEL',
    "instr_correction_input":'resp',
    "instr_correction_prefilt":[],
    "instr_correction_waterlevel":0.,
    
}

CONFIG_Preprocess = os.path.join('input','config_preprocess.json')

class ConfigPreprocess(object):
    """Contains basic parameters for the job (paths, etc.)"""

    def __init__(self):
        self.verbose = None
        self.testrun = None
        self.input_dirs = None
        self.input_format = None
        
        self.quality_minlengthsec = None
        self.quality_maxgapsec = None
        
        self.event_exclude = None
        self.event_exclude_winsec = None
        self.event_exclude_std = None
        self.event_exclude_n = None
        self.event_exclude_freqmin = None
        self.event_exclude_freqmax = None
        self.event_exclude_level = None
        
        self.wins = None
        self.wins_trim = None
        self.wins_detrend = None
        self.wins_demean = None
        self.wins_taper = None
        self.wins_cap = None
        self.wins_cap_threshold = None
        
        self.Fs_old = None
        self.Fs_new = None
        self.Fs_antialias_factor = None
        
        self.instr_correction = None
        self.instr_correction_unit = None
        self.instr_correction_input = None
        self.instr_correction_prefilt= None
        self.instr_correction_waterlevel = None
        
        
        self.initialize()
        
        
        
    def initialize(self):
        """Populates the class from ./config.json.
        If ./config.json does not exist, writes a default file and exits.
        """

        if not os.path.exists(CONFIG_Preprocess):
            with io.open(CONFIG_Preprocess, 'wb') as fh:
                json.dump(DEFAULT_Preprocess, fh, sort_keys=True, indent=4, separators=(",", ": "))
            return()

        # Load all options.
        with io.open(CONFIG_Preprocess, 'r') as fh:
            data = json.load(fh)
            
        for key, value in data.iteritems():
            setattr(self, key, value)
        
        # Make sure freqs. for downsampling are in descending order.
        self.Fs_new.sort() # Now in ascending order
        self.Fs_new=self.Fs_new[::-1] # Now in descending order


DEFAULT_Correlation = {
    "indirs": [],
    "time_begin": "2000-01-01T00:00:00.0000",
    "time_end": "2001-01-01T00:00:00.0000",
    "time_window_length":3600,
    "time_overlap":0,
    "time_min_window":3600,
    "corr_autocorr": False,
    "corr_type": "ccc",
    "corr_maxlag": 0,
    "corr_tensorcomponents": ["ZZ"],
    "n_stationpairs": 1,
    "input_format": "MSEED",
    "locations_mix": False,
    "locations":[],
    "whiten": False,
    "onebit": False,
    "ram_norm": False,
    "update": False,
}

CONFIG_Correlation = os.path.join('input','config_correlation.json')

class ConfigCorrelation(object):
    """Contains basic parameters for the correlation job (paths, etc.)"""

    def __init__(self):
        
        self.indirs = None
        
        self.time_begin = None
        self.time_end = None
        self.time_overlap = None
        self.time_window_length = None
        self.time_min_window = None
        self.corr_type = None
        self.corr_maxlag = None
        self.corr_tensorcomponents = None
        self.corr_autocorr = None
        self.input_format = None
        self.n_stationpairs = None
        self.locations_mix = None
        self.locations = None
        self.update = None
        self.whiten = None
        self.ram_norm = None
        self.onebit = None

        self.initialize()


    def initialize(self):

        """Populates the class from ./config.json.
        If ./config.json does not exist, writes a default file and exits."""

        if not os.path.exists(CONFIG_Correlation):
            with io.open(CONFIG_Correlation, 'wb') as fh:
                json.dump(DEFAULT_Correlation, fh, sort_keys=True, indent=4, separators=(",", ": "))
            return()

        # Load all options.
        with io.open(CONFIG_Correlation, 'r') as fh:
            data = json.load(fh)
            
        for key, value in data.iteritems():
            setattr(self, key, value)

