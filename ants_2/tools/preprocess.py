from __future__ import print_function

from numpy import isnan, isinf
from scipy.signal import sosfilt
from obspy import Stream, Inventory

# cfg should not be known here.
#from ants_2.config import ConfigPreprocess
#cfg = ConfigPreprocess()

class Preprocess(object):
    
    def __init__(self,stream,verbose,ofid):
        
        self.stream = stream
        if not isinstance(self.stream,Stream):
            msg = "stream must be an obspy stream object."
            raise TypeError(msg)
        self.verbose = verbose
        self.ofid = ofid
        
    
    
    def check_nan_inf(self):
        """
        Check if trace contains nan, inf and takes them out of the stream
        """
    
        for i in range(len(self.stream)):
            
            #- check NaN
            if True in isnan(self.stream[i].data):
                if self.verbose: 
                    print('** trace contains NaN, discarded',\
                file=self.ofid)
                
                del_trace = self.stream.pop(i)
                print(del_trace,file=self.ofid)
                continue
                        
            #- check infinity
            if True in isinf(self.stream[i].data):
                if self.verbose: print('** trace contains infinity, discarded',\
                file=self.ofid)
                
                del_trace = self.stream.pop(i)
                print(del_trace,file=self.ofid)
                continue
               

    def detrend(self):
        """
        remove linear trend
        """

        if self.verbose:
            print('* detrend\n',file=self.ofid)
        
        self.stream.detrend('linear')
    
    

    def demean(self):
        """
        remove the mean
        """

        if self.verbose: 
            print('* demean\n',file=self.ofid)
        
        self.stream.detrend('demean')

       
    
    def taper(self,perc):
        
        if self.verbose:
            print('* tapering\n',file=self.ofid)
            
        self.stream.taper(type='cosine',max_percentage=perc)
    
    def downsampling(self,Fs_new,sos_aa,zerophase_antialias=False):
        
        
        # Apply antialias filter
        
        for trace in self.stream:
        
            if zerophase_antialias:
                firstpass = sosfilt(sos_aa,trace.data)
                trace.data = sosfilt(sos_aa,firstpass[::-1])[::-1]
            else:
                trace.data = sosfilt(sos_aa,trace.data)
        
        # Decimate if possible, otherwise interpolate
        for i in range(len(Fs_new)):
            
            Fs_old = self.stream[0].stats.sampling_rate
            Fs = Fs_new[i]
            
            if Fs_old % Fs == 0:
                dec = int( Fs_old / Fs)
                self.stream.decimate(dec, no_filter=True, strict_length=False)
                if self.verbose:
                    print('* decimated traces to %g Hz' %Fs,
                    file=self.ofid)
            else:
                try:
                    self.stream.interpolate(sampling_rate = Fs,
                    method='lanczos')
                    print('* interpolated traces to %g Hz' %Fs,
                    file=self.ofid)
                except:
                    self.stream.interpolate(sampling_rate = Fs)
                    print('* interpolated trace to %g Hz' %Fs,
                    file=self.ofid)
       
            
   
    
    def remove_response(self,inv,pre_filt,**kwargs):
        
        """
        kwargs: Will be passed to obspy simulate or remove_response: water_level, taper,
        taper_fraction, plot, output
        """
        
        if isinstance(self.stream.inv,dict):
            self.stream.simulate(paz_remove=None,pre_filt=pre_filt,
            seedresp=inv,sacsim=True,pitsasim=False,**kwargs)
            if self.verbose:
                print('* removed instrument response using seedresp',file=self.ofid)
                
        elif isinstance(self.stream.inv,Inventory):
            self.stream.remove_response(inventory=inv,
            pre_filt=pre_filt,**kwargs)
            if self.verbose:
                print('* removed instrument response using staxml inv',file=self.ofid)
                
        else:
            msg = 'No inventory or seedresp found.'
            raise ValueError(msg)
        
        
    
    
