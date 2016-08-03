from __future__ import print_function

from numpy import isnan, isinf
from scipy.signal import sosfilt
from obspy import Trace, Inventory

# cfg should not be known here.
#from ants_2.config import ConfigPreprocess
#cfg = ConfigPreprocess()

class Preprocess(object):
    
    def __init__(self,trace,verbose,ofid):
        
        self.trace = trace
        if not isinstance(self.trace,Trace):
            msg = "trace must be an obspy trace object."
            raise TypeError(msg)
        self.verbose = verbose
        self.ofid = ofid
        
    
    
    def check_nan_inf(self):
        """
        Check if trace contains nan, inf
        """
        #- check NaN
        if True in isnan(self.trace.data):
            if self.verbose: print('** trace contains NaN, discarded',\
            file=self.ofid)
            return False
                    
        #- check infinity
        if True in isinf(self.trace.data):
            if self.verbose: print('** trace contains infinity, discarded',\
            file=self.ofid)
            return False
                    
        if self.verbose: 
            print('* number of points: '+str(self.
            trace.stats.npts)+\
        '\n',file=self.ofid)
        
        return True


    def detrend(self):
        """
        remove linear trend
        """

        if self.verbose:
            print('* detrend\n',file=self.ofid)
        
        self.trace.detrend('linear')
    
        return True


    def demean(self):
        """
        remove the mean
        """

        if self.verbose: 
            print('* demean\n',file=self.ofid)
        
        self.trace.detrend('demean')

        return True
    
    def downsampling(self,zerophase_antialias=False):
        
        
        # Apply antialias filter
        if zerophase_antialias:
            firstpass = sosfilt(self.trace.sos_aa,self.trace.data)
            self.trace.data = sosfilt(self.trace.sos_aa,firstpass[::-1])[::-1]
        else:
            self.trace.data = sosfilt(self.trace.sos_aa,self.trace.data)
        
        # Decimate if possible, otherwise interpolate
        for i in range(len(self.trace.Fs_new)):
            
            if self.trace.stats.Fs % self.trace.Fs_new[i] == 0:
                dec = int(self.trace.stats.Fs / self.trace.Fs_new[i])
                self.trace.decimate(dec, no_filter=True, strict_length=False)
                if self.verbose:
                    print('* decimated trace to %g Hz' %self.trace.Fs_new[i],
                    file=self.ofid)
            else:
                try:
                    self.trace.interpolate(sampling_rate = self.trace.Fs_new[i],
                    method='lanczos')
                    print('* interpolated trace to %g Hz' %self.trace.Fs_new[i],
                    file=self.ofid)
                except:
                    self.trace.interpolate(sampling_rate = self.trace.Fs_new[i])
                    print('* interpolated trace to %g Hz' %self.trace.Fs_new[i],
                    file=self.ofid)
        return True
            
   
    
    def remove_response(self,**kwargs):
        
        """
        kwargs: Will be passed to obspy simulate or remove_response: water_level, taper,
        taper_fraction, plot, output
        """
        
        if isinstance(self.trace.inv,dict):
            self.trace.simulate(paz_remove=None,pre_filt=self.trace.pre_filt,
            seedresp=self.trace.inv,sacsim=True,pitsasim=False,**kwargs)
            if self.verbose:
                print('* removed instrument response using seedresp',file=self.ofid)
                
        elif isinstance(self.trace.inv,Inventory):
            self.trace.remove_response(inventory=self.trace.inv,
            pre_filt=self.trace.pre_filt,**kwargs)
            if self.verbose:
                print('* removed instrument response using staxml inv',file=self.ofid)
                
        else:
            msg = 'No inventory or seedresp found.'
            raise ValueError(msg)
        
        return True
            