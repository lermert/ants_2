from __future__ import print_function
import io
import json
import os
import click


DEFAULT_Download = {
    "verbose":"true",
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
            print(data)
        for key, value in data.iteritems():
            setattr(self, key, value)

