#!/usr/bin/env python
"""
#####################################################################
# Author: Daniel Argueso <daniel> @ UIB
# Date:   2018-02-13T10:26:27+11:00
# Email:  d.argueso@uib.es
# Last modified by:   daniel
# Last modified time: 2018-02-13T10:26:29+11:00
#
# @Project@
# Version: x.0 (Beta)
# Description: Utilities to deal with WRF outputs
#
# Dependencies:
#
# Files:
#
#####################################################################
"""

import datetime as dt
import netCDF4 as nc
import numpy as np

def sel_wrfout_files(filelist,sdate,edate):
    """ Module to select files from a file list that have records between two
        given dates
        ----
        Output: list of files
    """

    d1 = dt.date(np.int(sdate[0:4]),np.int(sdate[5:7]),np.int(sdate[8:10]))
    d2 = dt.date(np.int(edate[0:4]),np.int(edate[5:7]),np.int(edate[8:10]))

    years=np.array([fname.split("/")[-1].split("_")[2][0:4] for fname in filelist],np.int)
    months=np.array([fname.split("/")[-1].split("_")[2][5:7] for fname in filelist],np.int)
    days=np.array([fname.split("/")[-1].split("_")[2][8:10] for fname in filelist],np.int)

    file_dates = np.array([dt.date(years[i],months[i],days[i]) for i in range(len(years))])

    selec_files=[filelist[i] for i,n in enumerate(file_dates) if  ((n>= d1) & (n<= d2))]

    return selec_files


###########################################################
###########################################################

def wrftime2date(files):
    """
    Conversion of dates from a wrf file or a list of wrf files
    format: [Y] [Y] [Y] [Y] '-' [M] [M] '-' [D] [D] '_' [H] [H] ':' [M] [M] ':' [S] [S]
    to a datetime object.
    """

    if len(files)==1:
        fin=nc.Dataset(str(files[0]),'r')
        times=fin.variables['Times']
    else:
        fin=nc.MFDataset(files[:])
        times=fin.variables['Times']

    year=np.zeros(len(times),dtype=np.int64)
    month=year.copy()
    day=year.copy()
    hour=year.copy()
    minute=year.copy()
    second=year.copy()

    for i in xrange(len(times)):
        listdate=times[i]
        year[i]=int(listdate[0])*1000 + int(listdate[1])*100 + int(listdate[2])*10 + int(listdate[3])
        month[i]=int(listdate[5])*10 + int(listdate[6])
        day[i]=int(listdate[8])*10 + int(listdate[9])
        hour[i]=int(listdate[11])*10 + int(listdate[12])
        minute[i]=int(listdate[14])*10 + int(listdate[15])
        second[i]=int(listdate[17])*10 + int(listdate[18])

    dates = [dt.datetime(year[i], month[i], day[i], hour[i], minute[i], second[i]) for i in xrange(len(times))]
    return dates
