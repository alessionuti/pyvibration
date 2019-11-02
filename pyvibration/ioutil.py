import re
import datetime

import numpy as np 
import pandas as pd 
import tkinter as tk
import tkinter.filedialog

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
        if re.search(r"([a-d])|([f-z])|([A-D])|([F-Z])", line): # ignore header lines
            continue
        else:
            headerlines = lineindex
            break
    dataf = pd.read_csv(file_path,r"\s+",header=None,skiprows=headerlines)
    infile.close()
    a = dataf.values
    return a


def writetxt_from_xlsx(xlsx_path,out_path,filenames):
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
        dataf = pd.read_excel(xlsx,sheet_name=sheetname)
        datarray = dataf.values
        columnlabels = list(dataf.columns.values)
        for i in range(1,len(columnlabels)):
            txt_array = np.column_stack((datarray[:,0],datarray[:,i]))
            headerstr = (sheetname+"\nCreated on "+str(datetime.datetime.utcnow().isoformat())
            +"\n"+columnlabels[0]+"\t"+columnlabels[i])
            outname = out_path+sheetname+"_"+filenames[i-1]+".txt"
            np.savetxt(outname,txt_array,fmt="%.10e", delimiter="\t",header=headerstr)


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
    filename = tkinter.filedialog.askopenfilename(parent=root,title=windowtitle)
    return filename


