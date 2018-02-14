#!/usr/bin/env python
""" Get_WRFvar Program that gets a variable from a list of wrfout and create a netcdf

Usage:
    get_wrfvar.py [-s SM] [-e EM] [-x WRUN] [-i HOST] [-p PATT] [-d DOM] VARNAME SY EY PATHIN PATHOUT
    get_wrfvar.py (-h | --help)
    get_wrfvar.py --version

Arguments:
    VARNAME     Requested variable name
    SY          First year to process
    EY          Last year to process
    PATHIN      Directory with input files
    PATHOUT     Directory to write out files


Options:
    -h --help               Show this screen
    --version               Show versions
    -s SM --smonth=SM       Optional first month    [default: 1].
    -e EM --emonth=EM       Optional last month     [default: 12].
    -x WRUN --exp=WRUN      Optional experiment named   [default: WRF_TEST].
    -i HOST --inst=HOST     Optional institution name for output info   [default: WRF].
    -p PATT --pattern=PATT  Optional starting pattern for wrfout [default: wrfout].
    -d DOM --domain=DOM    Optional domain to process [default: d01].


"""


import netCDF4 as nc
import numpy as np
from glob import glob
import datetime as dt
import calendar
import os
import atmopy as ap
from docopt import docopt

import atmopy.compute_vars as cvars
from atmopy.wrf_utils import wrftime2date,sel_wrfout_files


if __name__ == '__main__':
    args = docopt(__doc__, version='Get_WRFvar 1.0')

    print (args)


###########################################################
############# USER MODIF ##################################

wrun = args['--exp']

institution=args['--inst']

#path_in = "/home/dargueso/WRF_OUT/REHIPRE/%s/out" %(wrun)# Location of WRF output files
#path_out = "/home/dargueso/postprocessed/REHIPRE/%s/" %(wrun)#

path_in = args['PATHIN']
path_out = args['PATHOUT']

patt    = args['--pattern']                # Pattern of WRF output files



syear = int(args['SY'])
eyear = int(args['EY'])
smonth = int(args['--smonth'])
emonth = int(args['--emonth'])

varname = args['VARNAME']

dom = args['--domain']

###########################################################
###########################################################

y = syear
m = smonth
d = 1


if not glob('%s/%s_%s*' %(path_in,patt,dom)):
    raise Exception("ERROR: no available files in requested directory: %s/%s_%s*" %(path_in,patt,dom))

while (y < eyear or (y == eyear and m <= emonth)):
    print y, m, d

    #print "Extracting variables for day %s-%s-%s" %(y,str(m).rjust(2,"0"),str(d).rjust(2,"0"))

    if not glob('%s/%s_%s_%s*' %(path_in,patt,dom,y)):
        raise Exception("ERROR: no available files for year %s requested directory: %s/%s_%s_%s*" %(y,path_in,patt,dom,y))


    sdate="%s-%s-%s" %(y,str(m).rjust(2,"0"),str(d).rjust(2,"0"))


    filesin = sorted(glob('%s/%s_%s_%s*' %(path_in,patt,dom,sdate)))



    if not filesin:
        print "No available files for day %s" %(sdate)
        continue

    x=[]
    t=[]

    for n,filename in enumerate(filesin):
        print filename
        #print "Processing file %s of %s..." %(n, len(filesin))

        xFragment,atts = cvars.compute_WRFvar(filename,varname)
        if len(xFragment.shape) == 3:
            xFragment=np.expand_dims(xFragment,axis=0)
        if len(xFragment.shape) == 2:
            xFragment=np.expand_dims(xFragment,axis=0)
        tFragment = wrftime2date(filename.split())[:]

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

    edate = dt.datetime(y,m,d) + dt.timedelta(days=1)

    y = edate.year
    m = edate.month
    d = edate.day
