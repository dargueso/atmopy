#!/usr/bin/env python
"""
#####################################################################
# Author: Daniel Argueso <daniel> @ UIB
# Date:   2018-02-13T10:31:17+11:00
# Email:  d.argueso@uib.es
# Last modified by:   daniel
# Last modified time: 2018-02-13T10:32:07+11:00
#
# @Project@ atmopy
# Version: 0.1 (Beta)
# Description: Methods to calculate time-related variables.
#
# Dependencies: numpy, netCDF4, datetime
#
# Files:
#
#####################################################################
"""

import netCDF4 as nc
import numpy as np
import datetime as dt

def calc_nonstandard_seas(dates,nameseason):

  months_clim=np.array([dates[i].month for i,n in enumerate(dates)],np.int)
  monthsinit='JFMAMJJASONDJFMAMJJASOND'
  monthsseas=months_clim*0


  if nameseason=='Annual' or nameseason=='':
    monthsseas[:]=1

  else:
    startseas=monthsinit.index(nameseason)
    endseas=(monthsinit.index(nameseason)+len(nameseason))

    seaslen=[(i)%12+1 for i in range(startseas,endseas)]


    for x in range(len(seaslen)):
      monthsseas[months_clim==seaslen[x]]=1



  return monthsseas



def calc_seasons(dates):
    """ Function to calculate the season of each timestep provided a list of dates
    """
    months_clim=np.array([dates[i].month for i,n in enumerate(dates)],np.int)
    season_clim=months_clim.copy()
    season_clim[(months_clim==1) | (months_clim==2) | (months_clim==12)]=1
    season_clim[(months_clim==3) | (months_clim==4) | (months_clim==5)]=2
    season_clim[(months_clim==6) | (months_clim==7) | (months_clim==8)]=3
    season_clim[(months_clim==9) | (months_clim==10) | (months_clim==11)]=4

    return season_clim



def pptime2date(files):
    """
    Reading of dates from a file or a list of of postprocessed files and write
    to a datetime object.
    """
    if len(files)==1:
        fin=nc.Dataset(str(files[0]),'r')
        times=fin.variables['time']
        dates=nc.num2date(times[:],units=times.units,calendar=times.calendar)
    else:
        fin=nc.MFDataset(files)
        times=fin.variables['time']

        dates = nc.num2date(times[:],units=times.units,calendar=times.calendar)
    return dates


def is_month(month,smonth,emonth):
    if smonth>emonth:
        return (month >= smonth) | (month <= emonth)
    elif smonth<=emonth:
        return  (month >= smonth) & (month <= emonth)
