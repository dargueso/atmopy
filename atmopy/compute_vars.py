#!/usr/bin/env python
"""
#####################################################################
# Author: Daniel Argueso <daniel> @ UIB
# Date:   2017-11-28T18:25:06+01:00
# Email:  d.argueso@uib.es
# Last modified by:   daniel
# Last modified time: 2017-11-28T18:25:21+01:00
#
# @Project@ REHIPRE
# Version: x.0 (Beta)
# Description: Functions to calculate diagnostics from WRF outputs
#
# Dependencies: wrf-python, netCDF4, numpy
#
# Files:
#
#####################################################################
"""

import netCDF4 as nc
import numpy as np

def compute_PR(filename):
    """Function to calculate precipitation flux from a wrf output
       It also provides variable attribute CF-Standard
    """

    wrfvames = ['PREC_ACC_NC','PREC_ACC_C']
    varvals = {}
    ncfile = nc.Dataset(filename,'r')

    ## Specific to PR
    accum_dt = getattr(ncfile,'PREC_ACC_DT')

    ## Extracting variables required for diagnostic
    for var in wrfvames:
        varvals[var] = ncfile.variables[var][:]

    ## Computing diagnostic
    pr_acc = varvals['PREC_ACC_NC']+varvals['PREC_ACC_C']

    ## Deacumulating over prac_acc_dt (namelist entry)
    pr = pr_acc/(accum_dt*60.)



    atts = {"standard_name": "precipitation_flux",
                    "long_name"    : "total precipitation flux",
                    "units"        : "kg m-2 s-1"}

    ncfile.close()

    return pr,atts

def compute_PRACC(filename):
    """Function to calculate precipitation flux from a wrf output
       It also provides variable attribute CF-Standard
    """

    wrfvames = ['RAINNC','RAINC']
    varvals = {}
    ncfile = nc.Dataset(filename,'r')

    ## Extracting variables required for diagnostic
    for var in wrfvames:
        varvals[var] = ncfile.variables[var][:]

    ## Computing diagnostic
    pracc = varvals['RAINC']+varvals['RAINNC']

    atts = {"standard_name": "precipitation_accumulated",
                    "long_name"    : "total accumulated precipitation",
                    "units"        : "kg m-2"}

    ncfile.close()

    return pracc,atts

def compute_TAS(filename):
    """ Function to calculate 2-m temperature from WRF OUTPUTS
        It also provides variable attributes CF-Standard
    """

    ncfile = nc.Dataset(filename,'r')

    t2 = ncfile.variables['T2'][:]

    atts = {"standard_name": "air_temperature",
            "long_name":  "Surface air temperature",
            "units"    :  "K"                      ,
            "hgt"       :  "2 m"                    ,
            }

    return t2,atts



def compute_PSL(filename):
    """ Function to calculate PSL using wrf-python diagnostics
        It also provides variable attribute CF-Standard
        Note: this function was also coded manually for in compute_vars.py
    """

    from wrf import getvar,smooth2d

    ncfile = nc.Dataset(filename,'r')

    # Get the sea level pressure using wrf-python
    psl = getvar(ncfile, "slp")

    # Smooth the sea level pressure since it tends to be noisy near the mountains
    smooth_psl = smooth2d(psl, 3)

    atts = {"standard_name": "air_pressure_at_mean_sea_level",
                 "long_name":      "Sea Level Pressure"      ,
                 "units"    :      "Pa"                      ,
                }

    ncfile.close()
    return smooth_psl,atts


def compute_WA(filename):
    """ Function to calculate vertical wind W from WRF OUTPUTS
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    #Get vertical wind using wrf-python

    wa = getvar(ncfile, "wa")

    atts = {"standard_name": "vertical_wind_speed",
            "long_name":  "Vertical wind speed",
            "units"    :  "m s-1"                      ,
            "hgt"       :  ""                    ,
            }

    return wa,atts

def compute_WSFC(filename):
    """ Function to calculate vertical wind W from WRF OUTPUTS at the lowest level
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    #Get vertical wind using wrf-python

    wa = getvar(ncfile, "wa")

    if wa.ndim == 3:
        #No times
        wa = np.squeeze(wa[0,:,:])
    if wa.ndim == 4:
        #With times
        wa = np.squeeze(wa[:,0,:,:])

    atts = {"standard_name": "vertical_wind_speed",
            "long_name":  "Vertical wind speed",
            "units"    :  "m s-1"                      ,
            "hgt"       :  "lowest level"                    ,
            }

    return wa,atts


def compute_Wl2(filename):
    """ Function to calculate vertical wind W from WRF OUTPUTS at level 2
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    #Get vertical wind using wrf-python

    wa = getvar(ncfile, "wa")

    if wa.ndim == 3:
        #No times
        wa = np.squeeze(wa[1,:,:])
    if wa.ndim == 4:
        #With times
        wa = np.squeeze(wa[:,1,:,:])

    atts = {"standard_name": "vertical_wind_speed",
            "long_name":  "Vertical wind speed",
            "units"    :  "m s-1"                      ,
            "hgt"       :  "model level"                    ,
            }

    return wa,atts

def compute_uvmet10(filename):
    """ Function to calculate 10-m windspeed rotated to Earth coordinates
        from WRF OUTPUTS
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    #Get vertical wind using wrf-python

    wa = getvar(ncfile, "uvmet10")

    atts = {"standard_name": "wind_speed",
            "long_name":  "Earth coordinates wind speed",
            "units"    :  "m s-1"                      ,
            "hgt"       :  "10 m"                    ,
            }

    return wa,atts


def compute_uvmet10_wdir(filename):
    """ Function to calculate 10-m wind direction on Earth coordinates
        from WRF OUTPUTS
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    #Get wind direction (this function calculates both speed and direction)
    # We select direction only

    uvmet10_wind = getvar(ncfile, "uvmet10_wspd_wdir")
    atts = {"standard_name": "wind_direction",
            "long_name":  "Earth coordinates wind direction",
            "units"    :  "degrees_north"                      ,
            "hgt"       :  "10 m"                    ,
            }
    if uvmet10_wind.ndim == 3:
        #No times
        uvmet10_wind = np.squeeze(uvmet10_wind[1,:,:])


    return uvmet10_wind,atts

def compute_uvmet10_wspd(filename):
    """ Function to calculate 10-m wind speed on Earth coordinates
        from WRF OUTPUTS
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    #Get wind speed (this function calculates both speed and direction)
    #We select speed only

    uvmet10_wind = getvar(ncfile, "uvmet10_wspd_wdir")

    atts = {"standard_name": "wind_speed",
            "long_name":  "Earth coordinates wind speed",
            "units"    :  "m s-1"                      ,
            "hgt"       :  "10 m"                    ,
            }
    if uvmet10_wind.ndim == 3:
        #No times
        uvmet10_wind = np.squeeze(uvmet10_wind[0,:,:])


    return uvmet10_wind,atts

def compute_cloudfrac(filename):
    """ Function to calculate low med and high-leve cloud fraction
        from WRF OUTPUTS
        It also provides variable attributes CF-Standard
        Method obtained from NCL using pressure levels to determine the 3 categ.
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    cloudfrac = getvar(ncfile, "cloudfrac")

    atts = {"standard_name": "cloudfrac",
            "long_name":  "Cloud fraction",
            "units"    :  "%"                      ,
            "hgt"       :  ""                    ,
            }

    if cloudfrac.ndim == 3:
        #No times
        cloudfrac = np.squeeze(cloudfrac[0,:,:])

    return cloudfrac, atts

def compute_pressure(filename):
    """ Function to calculate pressure at Full model levels in hPa
        from WRF OUTPUTS.
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    pressure = getvar(ncfile, "pressure")

    atts = {"standard_name": "pressure",
            "long_name":  "full_model_level_pressure",
            "units"    :  "hPa"                      ,
            "hgt"       :  ""                    ,
            }

    return pressure, atts

def compute_height(filename):
    """ Function to calculate height at Full model levels in m
        from WRF OUTPUTS.
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    height = getvar(ncfile, "height")

    atts = {"standard_name": "height",
            "long_name":  "full_model_level_height",
            "units"    :  "m"                      ,
            "hgt"       :  ""                    ,
            }

    return height, atts
