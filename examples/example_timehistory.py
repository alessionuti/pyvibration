import numpy as np 
import pandas as pd 
import tkinter as tk
import matplotlib.pyplot as plt

import pyvibration
from pyvibration.ioutil import readtxtarray, openfiledialog

th = readtxtarray(openfiledialog("Select data file: time history"))

# Print some statistics
pyvibration.samplerate(th,printflag=True)
pyvibration.signal_stats(th)

# Plot the time history
plt.figure()
pyvibration.plots.thplot(th,"test","k")

# Plot the Fourier transform
plt.figure()
pyvibration.plots.fftplot(th,"test","b")

# Plot PSD estimates
plt.figure()
pyvibration.plots.psdplot_periodogram(th,"periodogram","#FFA500")
pyvibration.plots.psdplot_welch(th,"welch, 512 samples per segment","r",nperseg=512)
plt.ylim(1e-5,10)

plt.show()