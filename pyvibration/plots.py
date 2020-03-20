import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
from pyvibration import *


def psdplot(psd, name, lcolor, lsty="-", xmin=20, xmax=2000):
    """Plots PSD in log-log chart showing RMS
    """
    rmsval = rmslog(psd)
    legendstr = name + ", $G_{RMS}$=%3.2f" % rmsval
    plt.loglog(psd[:, 0], psd[:, 1], lsty, label=legendstr, color=lcolor)
    plt.legend()
    format_grid()
    plt.xlabel("Frequency [$Hz$]")
    plt.ylabel("PSD [$G^2/Hz$]")
    format_xaxis(plt.gca().xaxis, xmin, xmax)


def psdplot_from_frf(frf, psdin, name, lcolor, lsty="-", xmin=20, xmax=2000):
    """Plots PSD in log-log chart showing RMS, given FRF and input PSD
    """
    psd = frf2psd(frf, psdin, xmin, xmax)
    psdplot(psd, name, lcolor, lsty, xmin, xmax)


def psdplot_periodogram(sig, name, lcolor, lsty="-", xmin=20, xmax=2000):
    """Plots PSD (periodogram estimate) in log-log chart showing RMS, given a timeseries
    """
    psd = psd_periodogram(sig)
    psdplot(psd, name, lcolor, lsty, xmin, xmax)


def psdplot_welch(sig, name, lcolor, nperseg=None, lsty="-", xmin=20, xmax=2000):
    """Plots PSD (Welch estimate) in log-log chart showing RMS, given a timeseries
    """
    psd = psd_welch(sig, nperseg)
    psdplot(psd, name, lcolor, lsty, xmin, xmax)


def fftplot(sig, name, lcolor, decibel=False, lsty="-", xmin=20, xmax=2000):
    """Plots FFT, magnitude and phase, given a timeseries
    """
    legendstr = name
    afft = fft_real(sig)
    f = afft[:, 0]
    magnitude = np.abs(afft[:, 1] + 1j * afft[:, 2])
    phase = np.angle(afft[:, 1] + 1j * afft[:, 2]) * 180 / np.pi  # degrees

    plt.subplot(2, 1, 1)
    if (decibel):
        plt.semilogx(f, magnitude, lsty, label=legendstr, color=lcolor)
        plt.ylabel("Magnitude")
    else:
        plt.semilogx(f, 20 * np.log10(magnitude), lsty, label=legendstr, color=lcolor)
        plt.ylabel("Magnitude [$dB$]")
    plt.legend()
    format_grid()
    plt.xlabel("Frequency [$Hz$]")
    format_xaxis(plt.gca().xaxis, xmin, xmax)

    plt.subplot(2, 1, 2)
    plt.semilogx(f, phase, lsty, label=legendstr, color=lcolor)
    plt.legend()
    format_grid()
    plt.xlabel("Frequency [$Hz$]")
    plt.ylabel("Phase [$deg$]")
    plt.ylim(-180, 180)
    format_xaxis(plt.gca().xaxis, xmin, xmax)


def frfplot(frf, name, lcolor, lsty="-", xmin=20, xmax=2000):
    """Plots FRF magnitude in log-log chart
    """
    plt.loglog(frf[:, 0], frf[:, 1], lsty, label=name, color=lcolor)
    plt.legend()
    format_grid()
    plt.xlabel("Frequency [$Hz$]")
    plt.ylabel("|FRF|")
    format_xaxis(plt.gca().xaxis, xmin, xmax)


def thplot(sig, name, lcolor, lsty="-"):
    """Plots an acceleration time history
    """
    plt.plot(sig[:, 0], sig[:, 1], lsty, label=name, color=lcolor)
    plt.legend()
    format_grid()
    plt.xlabel("Time [$s$]")
    plt.ylabel("Acceleration [$G$]")


def format_xaxis(ax, xmin, xmax):
    """Format the log scale ticks on x axis
    """
    ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.ticklabel_format(style='plain', axis='x')
    plt.xticks(list(plt.xticks()[0]) + [xmin, xmax])
    plt.xlim(xmin, xmax)


def format_grid():
    """Format the grid
    """
    plt.grid(which="major", linestyle='-', linewidth=0.5)
    plt.grid(which="minor", linestyle='--', linewidth=0.5)
