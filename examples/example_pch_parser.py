import matplotlib.pyplot as plt

import pyvibration
from pyvibration.ioutil import read_pch_111

POINT_ID = 293979
SUBCASE = 1
COMPONENT = 3    # TZ

# Extract acceleration magnitude from punch file
frf = read_pch_111("./sample_punch.pch", SUBCASE, "ACCELERATION", POINT_ID, COMPONENT, "MAGNITUDE")

# Plot the frequency response
plt.figure()
pyvibration.plots.frfplot(frf, 'acceleration, POINT ID ' + str(POINT_ID), 'k')

plt.show()
