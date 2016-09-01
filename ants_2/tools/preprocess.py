import numpy as np
from obspy.signal.filter import envelope
from scipy import fftpack
from scipy.signal import iirfilter, zpk2sos


def bandpass(freqmin, freqmax, df, corners=4):
    """
    From obspy with modification.

    Butterworth-Bandpass Filter.

    Filter data from ``freqmin`` to ``freqmax`` using ``corners``
    corners.
    The filter uses :func:`scipy.signal.iirfilter` (for design)
    and :func:`scipy.signal.sosfilt` (for applying the filter).

    :type data: numpy.ndarray
    :param data: Data to filter.
    :param freqmin: Pass band low corner frequency.
    :param freqmax: Pass band high corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / order.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the filter order but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    low = freqmin / fe
    high = freqmax / fe
    # raise for some bad scenarios
    if high > 1:
        high = 1.0
        msg = "Selected high corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    if low > 1:
        msg = "Selected low corner frequency is above Nyquist."
        raise ValueError(msg)
    z, p, k = iirfilter(corners, [low, high], btype='band',
                        ftype='butter', output='zpk')
    sos = zpk2sos(z, p, k)
    return sos


def whiten_taper(freq1,freq2,df,npts,taper_width):
    
    #freqaxis=np.fft.fftfreq(tr.stats.npts,tr.stats.delta)
    
    ind_fw1 = int(round(freq1/df))
    ind_fw2 = int(round(freq2/df))
    
    
    length_taper = int(round((freq2-freq1)*\
    taper_width/df))
    
    taper_left = np.linspace(0.,np.pi/2,length_taper)
    taper_left = np.square(np.sin(taper_left))
    
    taper_right = np.linspace(np.pi/2,np.pi,length_taper)
    taper_right = np.square(np.sin(taper_right))
    
    taper = np.zeros(npts)
    taper[ind_fw1:ind_fw2] += 1.
    taper[ind_fw1:ind_fw1+length_taper] = taper_left
    taper[ind_fw2-length_taper:ind_fw2] = taper_right

    return taper


def whiten(tr,freq1,freq2,taper_width):
    # ToDo check fft here
    
    # Build a cosine taper for the frequency domain
    df = 1/(tr.stats.npts*tr.stats.delta)
    white_tape = whiten_taper(freq1,freq2,df,
        tr.stats.npts,taper_width)
    
    # Transform data to frequency domain
    tr.taper(max_percentage=0.05, type='cosine')
    spec = fftpack.fft(tr.data)
    
    # Don't divide by 0
    tol = np.max(np.abs(spec)) / 1e5
    
    # whiten
    spec /= np.abs(spec+tol)
    spec *= white_tape
    
    # Go back to time domain
    tr.data = np.real(fftpack.ifft(spec,n=len(tr.data)))
    #return tr
    
def cap(tr,cap_thresh):
    
    std = np.std(tr.data*1.e6)
    gllow = cap_thresh * std * -1
    glupp = cap_thresh * std
    tr.data = np.clip(tr.data*1.e6,gllow,glupp)/1.e6

    #return tr
    
def ram_norm(tr,winlen,prefilt=None):
    
    trace_orig = tr.copy()
    hlen = int(winlen*trace.stats.sampling_rate/2.)
    weighttrace = np.zeros(trace.stats.npts)
    
    if prefilt is not None:
        trace.filter('bandpass',freqmin=prefilt[0],freqmax=prefilt[1],\
        corners=prefilt[2],zerophase=True)
        
    envlp = envelope(trace.data)

    for n in xrange(hlen,trace.stats.npts-hlen):
        weighttrace[n] = np.sum(envlp[n-hlen:n+hlen+1]/(2.*hlen+1))
        
    weighttrace[0:hlen] = weighttrace[hlen]
    weighttrace[-hlen:] = weighttrace[-hlen-1]
    
    trace.data = trace_orig.data / weighttrace
    #return(trace_orig)