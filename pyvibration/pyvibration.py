import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.fftpack import fft
from scipy.signal import periodogram, welch


def interplog(x,y,xq,plotflag=False):
    """Logarithmic interpolation (base 10)
    Parameters
    ----------
    x : array
        Data x-values.
    y : array
        Data y-values.
    xq : array
        Interpolation x-values.
    plotflag : bool
        If True plots original and interpolated values. Default False.
    Returns
    -------
    yq : array
        Interpolated y-values.
    """
    logyq = np.interp(np.log10(xq),np.log10(x),np.log10(y))
    yq = 10**logyq
    if (plotflag):
        plt.figure()
        plt.loglog(x,y,'-r', label='raw')
        plt.loglog(xq,yq,'--k', label='interp')
        plt.title('raw vs log interp')
        plt.grid(which='both')
        plt.legend()
    return yq


def slopelog(a):  
    """Slopes of a loglog plot given the curve breakpoints
    Parameters
    ----------
    a : array
        Curve breakpoints table (x=a[:,0], y=a[:,1]).
    Returns
    -------
    slope : array
        Array of slopes.
    """ 
    x = a[:,0]
    y = a[:,1]

    n = len(x)      
    slope = np.zeros(n-1)
    for i in range(0,(n-1)):    
        slope[i] = np.log(y[i+1]/y[i])/np.log(x[i+1]/x[i])
    return slope


def rmslog(psd):
    """RMS value given the PSD curve breakpoints
    Parameters
    ----------
    psd : array
        PSD breakpoints (frequency=psd[:,0], PSD value=psd[:,1]).
    Returns
    -------
    rms : float
        RMS value.
    """ 
    x = psd[:,0]
    y = psd[:,1]
    
    n = len(x)
    slope = slopelog(psd)
    ra = 0
    for i in range (0,n-1):
        s = slope[i]
        if s < -1.0001 or s > -0.9999:
            ra += (y[i+1]*x[i+1]-y[i]*x[i])/(s+1)
        else:
            ra += y[i]*x[i]*np.log(x[i+1]/x[i])
    rms = np.sqrt(ra)
    return rms


def psd2frf(psdout,psdin,fmin,fmax,plotflag=0):
    """FRF amplitude given the PSD of input and response
    Parameters
    ----------
    psdout : array
        Response PSD breakpoints (frequency=psdout[:,0], PSD value=psdout[:,1]).
    psdin : array
        Input PSD breakpoints (frequency=psdin[:,0], PSD value=psdin[:,1]).
        Units shall be the same of psdout.
    fmin, fmax : float
        Minimum and maximum frequencies.
    plotflag : bool
        If True plots original and interpolated values for the interpolated array.
        Default False.
    Returns
    -------
    frf : array
        FRF amplitude breakpoints (frequency=frf[:,0], FRF amplitude=frf[:,1]).
    """
    freq, psdoutval, psdinval = interpolate_largest(psdout,psdin,plotflag)
    frfval = np.sqrt(psdoutval/psdinval)
    frf = truncate_x(np.column_stack((freq,frfval)), fmin, fmax)
    return frf


def frf2psd(frf,psdin,fmin,fmax,plotflag=0):
    """PSD of the response given the FRF amplitude and input PSD
    Parameters
    ----------
    frf : array
        FRF amplitude breakpoints (frequency=frf[:,0], FRF amplitude=frf[:,1]).
    psdin : array
        Input PSD breakpoints (frequency=psdin[:,0], PSD value=psdin[:,1]).
    fmin, fmax : float
        Minimum and maximum frequencies.
    plotflag : bool
        If True plots original and interpolated values for the interpolated array.
        Default False.
    Returns
    -------
    psdout : array
        Response PSD breakpoints (frequency=psdout[:,0], PSD value=psdout[:,1]).
        Units are the same of psdin.
    """
    freq, frfval, psdinval = interpolate_largest(frf,psdin,plotflag)
    psdoutval = ((frfval**2)*psdinval)
    psdout = truncate_x(np.column_stack((freq,psdoutval)), fmin, fmax)
    return psdout


def interpolate_largest(a1,a2,plotflag=0):
    """Interpolate two datasets on the longest x values array
    Parameters
    ----------
    a1, a2 : array
        Datasets (x=a[:,0], y=a[:,1]).
    plotflag : bool
        If True plots original and interpolated values for the interpolated array.
        Default False.
    Returns
    -------
    xout : array
        Common x-values.
    yout1 : array
        Fist set y-values.
    yout2 : array
        Second set y-values.
    """
    if (a1.shape[0] > a2.shape[0]):
        xout = a1[:,0]
        yout1 = a1[:,1]
        yout2 = interplog(a2[:,0],a2[:,1],xout,plotflag)
    else:
        xout = a2[:,0]
        yout1 = interplog(a1[:,0],a1[:,1],xout,plotflag)
        yout2 = a2[:,1]
    return xout, yout1, yout2


def truncate_x(a, xmin, xmax):
    """Truncates values with x external to a selected range
    Parameters
    ----------
    a : array
        Dataset (x=a[:,0], y=a[:,1]).
    xmin, xmax : float
        Minimum and maximum x-values.
    Returns
    -------
    aout : array
        Dataset (x=a[:,0], y=a[:,1]) with limited x range.
    """
    x = a[:,0]
    index = ((x > xmin) & (x < xmax))
    aout = a[index,:]
    return aout


def samplerate(sig,printflag=False):
    """Mean sample rate of a timeseries
    Parameters
    ----------
    sig : array
        Timeseries (t=sig[:,0], y=sig[:,1]).
    printflag: bool
        If True, prints min and max sample rate. Default False.
    Returns
    -------
    srmean : float
        Mean sample rate.
    Raises
    ------
    Exception
        for repeated time points or sample rate variation>1%
    """
    x = sig[:,0]
    
    num = len(x)
    dt = np.zeros(num-1)
    for i in range(0,num-1):
        dt[i]=x[i+1]-x[i]

    dtmin = dt.min()
    dtmax = dt.max()

    if(dtmin>1e-20):
        srmax = 1./dtmin
    else:
        srmax = np.Inf
        raise Exception("repeated time points")
    
    if(dtmax>1e-20):
        srmin = 1./dtmax
    else:
        srmin = np.Inf
        raise Exception("repeated time points")
    
    if(printflag):
        print ("dtmin = %8.4g s" % dtmin)
        print ("dtmax = %8.4g s \n" % dtmax)   
        print ("srmax = %8.4g Hz" % srmax)
        print ("srmin = %8.4g Hz \n" % srmin)

    srmean = np.mean(1./dt)

    if((srmax-srmin) > 0.01*srmean):
        raise Exception("variable sample rate (>1%)")
        
    return srmean


def fft_real(sig):
    """FFT of a real-valued timeseries
    Parameters
    ----------
    sig : array
        Timeseries (t=sig[:,0], y=sig[:,1]).
    Returns
    -------
    afft : array
        FFT (frequency=psd[:,0], FFT real part=psd[:,1], FFT imag part=psd[:,2]).
    """
    y = sig[:,1]
    sr = samplerate(sig)

    if ((len(y) % 2) != 0):     # make the length of timeseries even 
        y = y[0:-1]
    n = len(y)

    yfft = np.fft.rfft(y)       # fft
    yfft[1:-1] = 2*yfft[1:-1]   # double amplitude except zero & nyquist freq
    yfft /= n                   # normalization

    f = sr/n*np.arange(0,n/2+1) # frequency vector
    afft = np.column_stack((f,yfft.real,yfft.imag))
    return afft


def psd_periodogram(sig):
    """PSD estimate using periodogram
    Parameters
    ----------
    sig : array
        Timeseries (t=sig[:,0], y=sig[:,1]).
    Returns
    -------
    psd : array
        PSD estimate (frequency=psd[:,0], PSD value=psd[:,1]).
    """
    f,Pxx = periodogram(sig[:,1],samplerate(sig))
    psd = np.column_stack((f,Pxx))
    return psd


def psd_welch(sig,nperseg=None):
    """PSD estimate using Welch method with Hann windowing, 50% overlap
    Parameters
    ----------
    sig : array
        Timeseries (t=sig[:,0], y=sig[:,1]).
    nperseg: int
        Length of each segment. Defaults to 256.
    Returns
    -------
    psd : array
        PSD estimate (frequency=psd[:,0], PSD value=psd[:,1]).
    """
    f,Pxx = welch(sig[:,1],samplerate(sig),nperseg=nperseg)
    psd = np.column_stack((f,Pxx))
    return psd


def signal_stats(sig):
    """Print some statistical values of a timeseries
    Parameters
    ----------
    sig : array
        Timeseries (t=sig[:,0], y=sig[:,1]).
    """
    t = sig[:,0]
    y = sig[:,1]

    ave = np.mean(y)
    dur = t[-1]-t[0]
    rms = np.sqrt(np.mean(y**2))
    sd = np.std(y)
    skewness = stats.skew(y)
    kurtosis = stats.kurtosis(y,fisher=False)

    print ("          max = %8.4g " % max(y))
    print ("          min = %8.4g " % min(y))
    print ("         mean = %8.4g " % ave)
    print ("      std dev = %8.4g " % sd)
    print ("          rms = %8.4g " % rms)
    print ("     skewness = %8.4g " % skewness)
    print ("     kurtosis = %8.4g \n" % kurtosis)
    print ("        start = %8.4g s " % t[0])
    print ("          end = %8.4g s " % t[-1])
    print ("     duration = %8.4g s \n" % dur)


def psd_stats(psd):
    """Print RMS of acceleration, velocity, displacement from acceleration PSD
    Parameters
    ----------
    psd : array
        Acceleration PSD breakpoints (frequency=psd[:,0], PSD value=psd[:,1]).
    """
    f = psd[:,0]
    ya = psd[:,1]

    omega = 2*np.pi*f                   # rad/s
    
    yv = (ya*(9.8*1000)**2)/(omega**2)  # mm/s

    yd = yv/(omega**2)                  # mm
    
    rms = rmslog(psd)
    vrms = rmslog(np.column_stack((f,yv)))
    drms = rmslog(np.column_stack((f,yd)))
    
    print (" Acceleration ")
    print ("    RMS   = %8.4g g " % rms)
    print ("  3-sigma = %8.4g g \n" % (3*rms))
    print (" Velocity ") 
    print ("    RMS   = %8.4g mm/s " % vrms)
    print ("  3-sigma = %8.4g mm/s \n" % (3*vrms))
    print (" Displacement ") 
    print ("    RMS   = %8.4g mm " % drms)
    print ("  3-sigma = %8.4g mm \n" % (3*drms))
