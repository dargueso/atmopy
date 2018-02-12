#!/usr/bin/env python
"""wrf_time.py methods
Author: Daniel Argueso (d.argueso@unsw.edu.au) at CoECSS, UNSW, and CCRC, UNSW. Australia.

Methods to calculate time-related variables derived from WRF standard output. Standard NARCliM outputs.

Created: 08 August 2013


"""
import netCDF4 as nc
import numpy as np
import datetime as datetime

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

    dates = [datetime.datetime(year[i], month[i], day[i], hour[i], minute[i], second[i]) for i in xrange(len(times))]
    return dates

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
    