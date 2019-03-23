import matplotlib.cm as cm
import atmopy.obs_info as obs_info
import config as cfg
from atmopy.geo_info import get_res
import numpy as np
oinfo = obs_info.ObsInfo()

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans
from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm


def get_wrf_color_and_ticks(datasets,cmap='viridis',version_sep='HVC'):

    obs_datasets = oinfo.get_avail_datasets()
    res = list(set([get_res(wrun) for wrun in cfg.wrf_runs]))
    res.sort(reverse=True)
    c=cm.get_cmap(cmap)
    colors = list(iter(c(np.linspace(0,1,len(res)))))


    colorb=[]
    xtick =[]
    hatch =[]
    label =[]
    line={}
    obsn=0
    for dset in datasets:
        if dset in obs_datasets:
            colorb.append(np.asarray([0.7,0.7,0.7]))
            xtick.append(dset)
            hatch.append('')
            line[dset]='-'
            if obsn==0:
                label.append('OBS')
                obsn=1
            else:
                label.append('_nolegend_')

        else:
            colorb.append(colors[res.index(get_res(dset))])
            version = dset.split(version_sep)[-1]

            if version == '_SH':
                line[dset]='--'
                hatch.append('//')
                label.append('_nolegend_')
                xtick.append(str(get_res(dset)) + "km"+ "_SH")
            elif version == '_NC':
                line[dset]=':'
                hatch.append('/')
                label.append('_nolegend_')
                xtick.append(str(get_res(dset)) + "km"+ "_EX")
            else:
                line[dset]='-'
                hatch.append('')
                label.append(str(get_res(dset)) + "km")
                xtick.append(str(get_res(dset)) + "km"+ "_DP")

    return colorb,xtick,line,hatch,label


def add_map_mercator(lats,lons,truelat1=0.0):

    m = Basemap(projection='merc',llcrnrlat=lats[0],urcrnrlat=lats[1],llcrnrlon=lons[0],urcrnrlon=lons[1],lat_ts=truelat1,resolution='h',suppress_ticks=True,fix_aspect=True)
    m.drawcoastlines(linewidth=0.5, color='k')
    m.drawparallels(np.arange(-90.,90,5.0),labels=[1,0,0,0],linewidth=0.2,dashes=[3,3],fontsize=6)
    m.drawmeridians(np.arange(0.,360,5.0),labels=[0,0,0,1],linewidth=0.2,dashes=[3,3],fontsize=6)

    return m
def add_map_lambert(lats,lons,truelat1=30.,truelat2=90.,cen_lat=45.,cen_lon=0):

    m=Basemap(projection='lcc',llcrnrlat=lats[0],urcrnrlat=lats[1],llcrnrlon=lons[0],urcrnrlon=lons[1],\
    lat_1=truelat1,lat_2=truelat2,lat_0=cen_lat,lon_0=cen_lon,resolution='h', suppress_ticks=True,fix_aspect=True)
    m.drawcoastlines(linewidth=0.5, color='k')
    m.drawparallels(np.arange(-90.,90,5.0),labels=[1,0,0,0],linewidth=0.2,dashes=[3,3],fontsize=6)
    m.drawmeridians(np.arange(0.,360,5.0),labels=[0,0,0,1],linewidth=0.2,dashes=[3,3],fontsize=6)
    return m


def get_res_version(wrun,version_sep='HVC'):

    res = get_res(wrun)
    version = wrun.split(version_sep)[-1]

    return str(get_res(wrun)) + "km"+ version
