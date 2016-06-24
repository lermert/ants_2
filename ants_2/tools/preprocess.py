import numpy as np

from scipy.interpolate import interp1d
from obspy.core import Trace, Stream, UTCDateTime
from scipy.signal import cheby2,  cheb2ord,  filtfilt

from ants_2.config import ConfigPreprocess
cfg = ConfigPreprocess()

class Preprocess(object):
    
    def __init__(self,trace,ofid):
        
        self.trace = trace
        self.ofid = ofid
        
    
    
    def check_nan_inf(self):
    """
    Basic quality checks    
    """
        #- check NaN
        if True in np.isnan(self.trace.data):
            if cfg.verbose==True: print('** trace contains NaN, discarded',\
            file=self.ofid)
            return False
                    
        #- check infinity
        if True in np.isinf(trace.data):
            if cfg.verbose: print('** trace contains infinity, discarded',\
            file=self.ofid)
            return False
                    
        if cfg.verbose: print('* number of points: '+str(trace.stats.npts)+\
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
            print('* demean\n',file=ofid)
        
        self.trace.detrend('demean')

        return data


    def bandpass(self):
        """
        Butterworth bandpass
        """
        
        f_min = self.cfg.instr_correction_prefilt[0]
        f_max = self.cfg.instr_correction_prefilt[1]
        corners = self.cfg.instr_correction_prefilt[2]
        corners = self.cfg.instr_correction_prefilt[3]
        
        if self.verbose: 
            print('* bandpass between '+str(f_min)+' and '+str(f_max)+' Hz\n',
            file=self.ofid)
            print ('* filter order '+str(corners)+' \n',file=ofid)
        
        self.trace.filter('bandpass',freqmin=f_min,freqmax=f_max,
        corners=corners,zerophase=zerophase)
    
        return True

    
    
    def add_antialias(self):
        z, p, k = cheby2_lowpass(fs_old,freq)
        self.antialias = zpk2sos(z, p, k)
    
    def cheby2_lowpass(fs_old,freq)
        # From obspy
        nyquist = df * 0.5
        # rp - maximum ripple of passband, rs - attenuation of stopband
        rp, rs, order = 1, 96, 1e99
        ws = freq / nyquist  # stop band frequency
        wp = ws  # pass band frequency
        # raise for some bad scenarios
        if ws > 1:
            ws = 1.0
            msg = "Selected corner frequency is above Nyquist. " + \
                  "Setting Nyquist as high corner."
            warnings.warn(msg)
        while True:
            if order <= maxorder:
                break
            wp = wp * 0.99
            order, wn = cheb2ord(wp, ws, rp, rs, analog=0)
        return cheby2(order, rs, wn, btype='low', analog=0, output='zpk')
    
    
        
    def lowpass(self):
        
        # determine the lowest frequency, and from that the antialias filter frequency
        # filter
        # decimate in steps.
        
        if self.verbose:
            print('* lowpass below '+str(f_max)+' Hz\n',file=ofid)
        
        data.filter('lowpass', freq=f_max,corners=corners,zerophase=False)
    

        return data

    #==================================================================================================
    # REMOVE INSTRUMENT RESPONSE
    #==================================================================================================

    def remove_response(data,respdir,unit,freqs,waterlevel,verbose, ofid):

        """
        Remove instrument response located in respdir from data. Unit is displacement (DIS), velocity (VEL) or acceleration (ACC).

        Return 1 if successful. Return 0 otherwise.
        """

        #- RESP file ==================================================================================

        resp_file=respdir+'/RESP.'+data.stats.network+'.'+data.stats.station+'.'+data.stats.location+'.'+data.stats.channel

        if verbose==True:
            print('* RESP file: '+resp_file+'\n',file=ofid)

        #- try to remove response if the RESP file exists =============================================

        if os.path.exists(resp_file):

            success=1

            if verbose==True:
                print('* remove instrument response, unit='+unit+'\n',file=ofid)
                    
            resp_dict = {"filename": resp_file, "units": unit, "date": data.stats.starttime}
        
            try:
                data.simulate(seedresp=resp_dict, water_level=float(waterlevel),\
                nfft_pow2=True, simulate_sensitivity=False,pre_filt=tuple(freqs),\
                pitsasim=False,sacsim=True)
            except ValueError:
                if verbose==True: 
                    print('** could not remove instrument response\n',file=ofid)
                success=0

        #- response cannot be removed because RESP file does not exist

        else:
            if verbose==True: 
                print('** could not find correct RESP file\n',file=ofid)
            
            success=0

        return success, data


    
        
        

def preprocess_segment(cfg,trace,ofid):
                
    
                #===============================================================================
                # processing (detrending, filtering, response removal, decimation)
                #===============================================================================

    if cfg.detrend:
    
            trace=proc.detrend(trace,cfg.verbose,ofid)
            
                if check:
                    dfile.write('Detrended\n')
                    print(trace.data[0:20],file=dfile)
                    dfile.write('\n')
            
            if cfg.demean:
    
            trace=proc.demean(trace,cfg.verbose,ofid)
            
                if check:
                    dfile.write('Mean removed\n')
                    print(trace.data[0:20],file=dfile)
                    dfile.write('\n')
            
                if cfg.cap_glitches:
                    
                        std = np.std(trace.data/1.e6)
                    gllow = cfg.cap_threshold * -std
                    glupp = cfg.cap_threshold * std
                    trace.data = np.clip(trace.data/1.e6,gllow,glupp)*1.e6
                
         
         
    #- event exclusion based on energy levels.. ========================================================================                    
                # This should operate directly on the trace.
                if cfg.exclude_events:
                    ee.event_exclude(trace,cfg.exclude_windows,cfg.exclude_n,\
                    cfg.exclude_freq,cfg.exclude_level)
                
    
                #- taper edges ========================================================================
    
                if cfg.taper_do == True:
    
                    trace=proc.taper(trace,cfg.taper_width,cfg.verbose,ofid)
                
                    if check == True:
                        dfile.write('Tapered\n')
                        print(trace.data[0:20],file=dfile)
                        dfile.write('\n')
            
                #- downsampling =======================================================================
                sampling_rate_index=0
                while sampling_rate_index<len(Fs_new):
                    if trace.stats.sampling_rate>Fs_new[sampling_rate_index]:
                        trace=proc.downsample(trace,Fs_new[sampling_rate_index],\
                        cfg.verbose,ofid)
                    sampling_rate_index+=1
                newtrace = trace.copy()
                del trace
               
                if check == True:
                    dfile.write('(Downsampled), copied\n')
                    print(newtrace.data[0:20],file=dfile)
                    dfile.write('\n')   
                #- remove instrument response =========================================================
    
                if cfg.remove_response == True:
    
                    removed,newtrace=proc.remove_response(newtrace,respdir,unit,\
                    freqs,wl,cfg.verbose,ofid)
                    if removed==False:
                        print('** Instrument response could not be removed! \
                        Trace discarded.',file=ofid)
                        continue
                    
                    if True in np.isnan(newtrace):
                        print('** Deconvolution seems unstable! Trace discarded.',\
                        file=ofid)
                        continue

