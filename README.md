# pyvibration

A Python vibration data analysis library

## Description

The aim is to provide a convenient Python interface for vibration data analysis (fourier analysis, PSD estimation, plotting) using where possible numpy/scipy implementations, with the addition of some utilities.  

Data sets are represented as numpy arrays where the first column is time or frequency. 

Module ```ioutil``` contains functions to load data from text and excel files.

Module ```plots``` is used for plotting time histories and spectra.

## Prerequisites

* [NumPy](https://numpy.org/)
* [SciPy](https://www.scipy.org/)
* [Matplotlib](https://matplotlib.org/)
* [pandas](https://pandas.pydata.org/)
* package [nastran_pch_reader](https://pypi.org/project/nastran_pch_reader/ for parsing NASTRAN punch files

## Use

```
import pyvibration
```

See examples folder for some examples (did you ever imagine that?).

## License

This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

* Tom Irvine http://www.vibrationdata.com/, for the many useful notes and scripts about vibration and shock analysis.
* Nikolay Asmolovskiy https://github.com/anick107/, for the python pch parser.
