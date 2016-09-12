
import numpy as np
from math import sqrt, isnan
from obspy.signal.cross_correlation import xcorr
from scipy.signal.signaltools import fftconvolve
import warnings

def my_centered(arr, newsize):
    # get the center portion of a 1-dimensional array
    n = len(arr)
    i0 = (n - newsize) // 2
    if n%2 == 0:
        i0 += 1
    i1 = i0 + newsize
    return arr[i0:i1]

def classic_xcorr(trace1, trace2, max_lag_samples):
   
    x_corr = xcorr(trace1.data, trace2.data,\
        max_lag_samples, True)[2]
    
    return x_corr

def get_correlation_params(data1,data2):

    if len(data1) == 0 or len(data2) == 0:
        return(0,0,0,0,0,0)
    # Get the signal energy; most people normalize by the square root of that
    ren1 = np.correlate(data1,data1,mode='valid')[0]
    ren2 = np.correlate(data2,data2,mode='valid')[0]

    # Get the window rms
    
    rms1 = sqrt(ren1 / len(data1))
    
    rms2 = sqrt(ren2 / len(data2)) 
    
    
    # A further parameter to 'see' impulsive events: range of standard deviations
    nsmp = int(len(data1)/4)

    std1 = [0,0,0,0]
    std2 = [0,0,0,0]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in range(4):
    
            
            std1[i] = np.std(data1[i*nsmp:(i+1)*nsmp])
            if isnan(std1[i]):
                return(0,0,0,0,0,0)
            std2[i] = np.std(data2[i*nsmp:(i+1)*nsmp])
            if isnan(std1[i]):
                return(0,0,0,0,0,0)
   
    # Add small value not to divide by zero
    tol = np.max(std1) * 1e-6 
    if tol != 0:
        rng1 = max(std1) / (min(std1) + tol)
        rng2 = max(std2) / (min(std2) + tol)
    else:
        rng1 = 0
        rng2 = 0

    return(rms1,rms2,ren1,ren2,rng1,rng2)

    
def cross_covar(data1, data2, max_lag_samples, normalize, params=False):
    
    #ToDo: deal with params
    
# remove mean and normalize; this should have no effect on the energy-normalized 
#correlation result, but may avoid precision issues if trace values are very small
    #if normalize:
    #    scale1 = 1./np.max(np.abs(data1))
    #    scale2 = 1./np.max(np.abs(data2))
    #    data1*=scale1
    #    data2*=scale2
    
    if len(data1) == 0 or len(data2) == 0:
        return([],[])

    
    data1-=np.mean(data1)
    data2-=np.mean(data2)
        
    # Make the data more convenient for C function np.correlate

    data1 = np.ascontiguousarray(data1, np.float32)
    data2 = np.ascontiguousarray(data2, np.float32)
    
    if params:
        params = get_correlation_params(data1,data2)
        ren1, ren2 = params[2:4]
    else:
        ren1 = np.correlate(data1,data1,mode='valid')[0]
        ren2 = np.correlate(data2,data2,mode='valid')[0]

    if ren1 == 0.0 or ren2 == 0.0 and normalize:
        return([],[])



    # scipy.fftconvolve is way faster than np.correlate, and zeropads for non-circular convolution
    ccv = fftconvolve(data1[::-1],data2,mode='same')

    #if normalize:
    #    ccv /= (scale1*scale2) 
    


    if normalize:
        ccv /= ( sqrt(ren1) * sqrt(ren2) )

    return my_centered(ccv,2*max_lag_samples+1),params



