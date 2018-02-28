#!/usr/bin/env python
""" Get_WRFvar Program that gets a variable from a list of wrfout and create a netcdf

"""


import netCDF4 as nc
import numpy as np
from glob import glob
import datetime as dt
import calendar
import os
import atmopy as ap
from optparse import OptionParser
import re

import atmopy.compute_vars as cvars
from atmopy.wrf_utils import wrftime2date,sel_wrfout_files


def read_input(filename):
    """
    Read input file with input arguments:

    wrun: Simulation ID - used in output path and naming
    institution: Insititution - used in output files
    path_in: path to input wrf files
    path_out: path to output postprocessed files

    patt: pattern of input files - generally wrfout or wrfhrly

    syear: First year to postprocess (e.g. 1990)
    eyear: Last year to postprocess (e.g. 2009)
    smonth: First month to postprocess (e.g. 1)
    emonth: Last month to postprocess (e.g. 12)

    varname: Variable to postprocess (e.g. PR)

    dom: domain to postprocess (e.g. 'd04')

    acc_dt: accumuluated period for variables such as precipitation (in minutes)
    """

    filein=open(filename,'r')
    lines=filein.readlines()

    inputinfs={}
    entryname=[]
    entryvalue=[]
    for line in lines:
      line=re.sub('\s+',' ',line)
      line=line.replace(" ", "")
      li=line.strip()
      #Ignore empty lines
      if li:
        #Ignore commented lines
        if not li.startswith("#"):
          values=li.split('=')
          entryname.append(values[0])
          entryvalue.append(values[1])
    filein.close()

    for ii in xrange(len(entryname)):
      inputinfs[entryname[ii]]=entryvalue[ii]

    return inputinfs




#### READING INPUT FILE ######
### Options

parser = OptionParser()

parser.add_option("-i", "--infile", dest="infile",
help="file with the input arguments", metavar="INPUTFILE")
(opts, args) = parser.parse_args()

###


#### Reading input info file ######
inputinfs=read_input(opts.infile)


###########################################################
###########################################################

wrun = inputinfs['wrun']
institution=inputinfs['institution']
path_in = inputinfs['path_in']
path_out = inputinfs['path_out']
patt    = inputinfs['patt']
syear = int(inputinfs['syear'])
eyear = int(inputinfs['eyear'])
smonth = int(inputinfs['smonth'])
emonth = int(inputinfs['emonth'])
varname = inputinfs['varname']
dom = inputinfs['dom']

###########################################################
###########################################################

y = syear
m = smonth
d = 1


if not glob('%s/%s_%s*' %(path_in,patt,dom)):
    raise Exception("ERROR: no available files in requested directory: %s/%s_%s*" %(path_in,patt,dom))

while (y < eyear or (y == eyear and m <= emonth)):
    print y, m, d

    if not glob('%s/%s_%s_%s*' %(path_in,patt,dom,y)):
        raise Exception("ERROR: no available files for year %s requested directory: %s/%s_%s_%s-%s*" %(y,path_in,patt,dom,y,m))


    sdate="%s-%s" %(y,str(m).rjust(2,"0"))


    filesin = sorted(glob('%s/%s_%s_%s*' %(path_in,patt,dom,sdate)))



    if not filesin:
        print "No available files for month %s" %(sdate)
        continue

    x=[]
    t=[]

    if len(filesin) == 1:
        print filesin
        varout = cvars.compute_WRFvar(filesin[0],varname,inputinfs)
        otimes =  wrftime2date(filesin[0].split())[:]

    else:

        for n,filename in enumerate(filesin):
            print filename

            tFragment = wrftime2date(filename.split())[:]
            xFragment,atts = cvars.compute_WRFvar(filename,varname,inputinfs)

            if len(tFragment)==1:
                if len(xFragment.shape) == 3:
                    xFragment=np.expand_dims(xFragment,axis=0)
                if len(xFragment.shape) == 2:
                    xFragment=np.expand_dims(xFragment,axis=0)


            x.append(xFragment)
            t.append(tFragment)

        varout = np.concatenate(x, axis=0)
        otimes = np.concatenate(t, axis=0)



    ###########################################################
    ###########################################################

    # ## Creating netcdf files
    fileout = "%s/%s_%s_%s.nc" %(path_out,institution,varname,str(sdate))

    ref_file = nc.Dataset(filesin[0])
    lat=ref_file.variables['XLAT'][0,:]
    lon=ref_file.variables['XLONG'][0,:]

    varinfo = { 'values': varout,
                'varname': varname,
                'atts':atts,
                'lat': lat,
                'lon': lon,
                'times': otimes}

    cvars.create_netcdf(varinfo,fileout)


    #edate = dt.datetime(y,m,d) + dt.timedelta(days=1)
    print otimes[-1].strftime("%Y-%m-%d")
    edate = otimes[-1] + dt.timedelta(days=1)

    y = edate.year
    m = edate.month
    d = edate.day
