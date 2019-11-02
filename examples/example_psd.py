import numpy as np 
import pandas as pd 
import tkinter as tk
import matplotlib.pyplot as plt

import pyvibration
from pyvibration.ioutil import readtxtarray, openfiledialog

psdin = readtxtarray(openfiledialog("Select data file: input PSD"))
frf = readtxtarray(openfiledialog("Select data file: FRF amplitude"))

# Calculate response PSD
psdout = pyvibration.frf2psd(frf,psdin,20,2000)

# Print some statistics of the response
pyvibration.psd_stats(psdout)

# Plot the response PSD
plt.figure()
pyvibration.plots.psdplot(psdout,"response","b")

plt.show()