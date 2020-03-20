import re
import datetime

import numpy as np
import pandas as pd
import tkinter as tk
import tkinter.filedialog
import nastran_pch_reader


def readtxtarray(file_path):
    """Reads a text file automatically detecting header lines
    Parameters
    ----------
    file_path : string
        Input text file path.
    Returns
    -------
    a : array
        Array of the numeric data.
    """
    infile = open(file_path, encoding="utf-8", errors="replace")
    lines = infile.readlines()
    for lineindex, line in enumerate(lines):
        if re.search(r"([a-d])|([f-z])|([A-D])|([F-Z])", line):  # ignore header lines
            continue
        else:
            headerlines = lineindex
            break
    dataf = pd.read_csv(file_path, r"\s+", header=None, skiprows=headerlines)
    infile.close()
    a = dataf.values
    return a


def writetxt_from_xlsx(xlsx_path, out_path, filenames):
    """Reads data from an excel file and writes to multiple txt files
    The first column of each sheet must be the independent variable (time or frequency);
    the following columns are an arbitrary number of data columns.
    The output files are tab-separated text files (one for each data column) with two columns
    (time/freq, data); the header includes column names from excel.
    Output file names are [sheetname]_[filename].txt
    Parameters
    ----------
    xlsx_path: string
        Excel file path.
    out_path: string
        Output directory.
    filenames: list of string
        Output names list.
    """
    xlsx = pd.ExcelFile(xlsx_path)
    sheetnames = xlsx.sheet_names
    for sheetname in sheetnames:
        dataf = pd.read_excel(xlsx, sheet_name=sheetname)
        datarray = dataf.values
        columnlabels = list(dataf.columns.values)
        for i in range(1, len(columnlabels)):
            txt_array = np.column_stack((datarray[:, 0], datarray[:, i]))
            headerstr = (sheetname + "\nCreated on " + str(datetime.datetime.utcnow().isoformat())
                         + "\n" + columnlabels[0] + "\t" + columnlabels[i])
            outname = out_path + sheetname + "/" + filenames[i - 1] + ".txt"
            np.savetxt(outname, txt_array, fmt="%.10e", delimiter="\t", header=headerstr)


def openfiledialog(windowtitle):
    """Shows a dialog to select a file
    Parameters
    ----------
    windowtitle: string
        Dialog title.
    Returns
    -------
    filename : string
        Selected file path.
    """
    root = tk.Tk()
    root.withdraw()
    filename = tkinter.filedialog.askopenfilename(parent=root, title=windowtitle)
    return filename


def read_pch_111(file_path, subcase, request, eid, index, complex_type):
    """Reads a NASTRAN SOL 111 .pch file
    Parameters
    ----------
    file_path : string
        Input pch file path.
    subcase : int
        SOL 111 subcase number
    request : string
        allowed values:
            'ACCELERATION'
            'DISPLACEMENTS'
            'MPCF'
            'SPCF'
            'ELEMENT FORCES' only for CELAS2 and CBUSH elements
    eid : int
        POINT ID or ELEMENT ID based on request type
    index : int
        index of the component to return (1-based)
    complex_type : string
        allowed values:
            'MAGNITUDE'
            'PHASE'
            'REAL'
            'IMAGINARY'
    Returns
    -------
    a : array
        Array of the numeric data (frequency=a[:,0], FRF=a[:,1]).
    """
    # parse punch file and store into PchParser.parsed_data
    parser = nastran_pch_reader.PchParser(file_path)
    # extract variables and convert to numpy arrays
    freq = np.array(parser.get_frequencies(subcase))
    if request == 'ACCELERATION':
        values = np.array(parser.get_accelerations(subcase)[eid])
    elif request == 'DISPLACEMENTS':
        values = np.array(parser.get_displacements(subcase)[eid])
    elif request == 'MPCF':
        values = np.array(parser.get_mpcf(subcase)[eid])
    elif request == 'SPCF':
        values = np.array(parser.get_spcf(subcase)[eid])
    elif request == 'ELEMENT FORCES':
        values = np.array(parser.get_forces(subcase)[eid])
    else:
        raise KeyError('invalid request type')
    # extract component and compute output
    if complex_type == 'MAGNITUDE':
        val = np.abs(values[:, index - 1])
    elif complex_type == 'PHASE':
        val = np.angle(values[:, index - 1], deg=True)
    elif complex_type == 'REAL':
        val = np.real(values[:, index - 1])
    elif complex_type == 'IMAGINARY':
        val = np.imag(values[:, index - 1])
    else:
        raise KeyError('invalid complex type')
    return np.column_stack((freq, val))
