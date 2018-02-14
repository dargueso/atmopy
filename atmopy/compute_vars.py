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
import datetime as dt


def compute_WRFvar (filename,varname):
    """ Function to compute a variable name from wrf outputs
        filename: wrfout (or other file name used as input)
        varname : variable to be extracted or computed from WRF outputs
    """

    ncfile = nc.Dataset(filename,'r')



    if varname in ncfile.variables.keys():
        varval = ncfile.variables[varname][:]
        varatt={}
        for att in ncfile.variables[varname].ncattrs():
            varatt[att] = getattr(ncfile.variables[varname],att)

    else:
        #Old version when is out of the functtion
        #compute=getattr('compute_'+varname)
        method_name='compute_%s' %(varname)
        possibles = globals().copy()
        possibles.update(locals())
        compute=possibles.get(method_name)

        varval, varatt=compute(filename)

    ncfile.close()

    return varval,varatt



###########################################################
###########################################################

def create_netcdf(var,filename):
    print('\n Create output file %s') %(filename)

    otimes = var['times']
    outfile = nc.Dataset(filename,'w',format='NETCDF3_CLASSIC')

    outfile.createDimension('time',None)
    outfile.createDimension('bnds',2)
    if var['values'].ndim == 4:

        outfile.createDimension('y',var['values'].shape[2])
        outfile.createDimension('x',var['values'].shape[3])
        outfile.createDimension('lev',var['values'].shape[1])

        outvar  = outfile.createVariable(var['varname'],'f8',('time','lev','y','x'),fill_value=1e20)

    if var['values'].ndim == 3:
        outfile.createDimension('y',var['values'].shape[1])
        outfile.createDimension('x',var['values'].shape[2])

        outvar  = outfile.createVariable(var['varname'],'f8',('time','y','x'),fill_value=1e20)

    outtime = outfile.createVariable('time','f8','time',fill_value=1e20)
    outtime_bnds = outfile.createVariable('time_bnds','f8',('time','bnds'),fill_value=1e20)
    outlat  = outfile.createVariable('lat','f8',('y','x'),fill_value=1e20)
    outlon  = outfile.createVariable('lon','f8',('y','x'),fill_value=1e20)

    if var['values'].ndim == 4:
        outlev = outfile.createVariable('levels','f8',('lev'),fill_value=1e20)
        if var['varname']=='cloudfrac':
            setattr(outlev,"standard_name","cloud-level")
            setattr(outlev,"long_name","Clouds level")
            setattr(outlev,"units","")
            setattr(outlev,"_CoordinateAxisType","z")
        else:
            setattr(outlev,"standard_name","model-level")
            setattr(outlev,"long_name","Model level")
            setattr(outlev,"units","eta levels")
            setattr(outlev,"_CoordinateAxisType","z")

    setattr(outlat,"standard_name","latitude")
    setattr(outlat,"long_name","Latitude")
    setattr(outlat,"units","degrees_north")
    setattr(outlat,"_CoordinateAxisType","Lat")


    setattr(outlon,"standard_name","longitude")
    setattr(outlon,"long_name","Longitude")
    setattr(outlon,"units","degrees_east")
    setattr(outlon,"_CoordinateAxisType","Lon")

    setattr(outtime,"standard_name","time")
    setattr(outtime,"long_name","Time")
    setattr(outtime,"units","hours since 1949-12-01 00:00:00")
    setattr(outtime,"calendar","standard")

    setattr(outtime_bnds,"standard_name","time_bnds")
    setattr(outtime_bnds,"long_name","time_bounds")
    setattr(outtime_bnds,"units","hours since 1949-12-01 00:00:00")
    setattr(outtime_bnds,"calendar","standard")


    step_seconds = np.int((otimes[1]-otimes[0]).total_seconds())

    outtime[:] = nc.date2num([otimes[x] for x in range(len(otimes))],units='hours since 1949-12-01 00:00:00',calendar='standard')

    outtime_bnds[:,0]=nc.date2num([otimes[x]-dt.timedelta(seconds=step_seconds/2.) for x in range(len(otimes))],units='hours since 1949-12-01 00:00:00',calendar='standard')
    outtime_bnds[:,1]=nc.date2num([otimes[x]+dt.timedelta(seconds=step_seconds/2.) for x in range(len(otimes))],units='hours since 1949-12-01 00:00:00',calendar='standard')


    outlat[:]  = var['lat'][:]
    outlon[:]  = var['lon'][:]

    outvar[:] = var['values'][:]

    for outatt in var['atts'].keys():
        setattr(outvar,outatt,var['atts'][outatt])

    setattr(outfile,"creation_date",dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(outfile,'author','Daniel Argueso @UIB')
    setattr(outfile,'contact','d.argueso@uib.es')
    #setattr(outfile,'comments','files created from wrf outputs %s/%s' %(path_in,patt))

    outfile.close()

###########################################################
###########################################################

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
            "hgt"       :  "low_mid_high",
            "low_thres": "300 m agl",
            "mid_thres": "2000 m agl",
            "high_thres": "6000 m agl",

            }

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

def compute_TA(filename):
    """ Function to calculate temperature in kelvin from WRF outputs
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    tk = getvar(ncfile,"tk")

    atts = {"standard_name": "Temperature",
            "long_name":  "air_temperature",
            "units"    :  "K"                      ,
            "hgt"       :  ""                    ,
            }
    return tk,atts

def compute_UA(filename):
    """ Function to calculate earth-rotated Eastward wind components in m s-1 from WRF outputs
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    ua = np.squeeze(getvar(ncfile,"uvmet")[0,:])

    atts = {"standard_name": "eastward_wind",
            "long_name":  "Eastward Wind",
            "units"    :  "m s-1"                      ,
            "hgt"       :  ""                    ,
            }
    return ua,atts

def compute_VA(filename):
    """ Function to calculate earth-rotated Northward wind components in m s-1 from WRF outputs
        It also provides variable attributes CF-Standard
    """
    from wrf import getvar

    ncfile = nc.Dataset(filename,'r')

    va = np.squeeze(getvar(ncfile,"uvmet")[1,:])

    atts = {"standard_name": "northward_wind",
            "long_name":  "Northward Wind",
            "units"    :  "m s-1"                      ,
            "hgt"       :  ""                    ,
            }
    return va,atts
