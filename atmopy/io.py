
import xarray as xr
from atmopy.utils import is_month
import atmopy.obs_info as obs_info
from glob import glob
import pandas as pd
import numpy as np
import netCDF4 as nc

oinfo = obs_info.ObsInfo()



def load_wrf(wrun,freq,var,levs=[]):

    if levs != []:
        patt_in = cfg.patt_in + "_" + levs

    else:
        patt_in = cfg.patt_in

    filesin = sorted(glob('%s/%s/20??-20??/%s_%s_%s_*.nc' %(cfg.path_in,wrun,patt_in,freq,var)))
    #print ('%s/%s/20??-20??/%s_%s_%s_*.nc' %(cfg.path_in,wrun,patt_in,freq,var))
    fin = xr.open_mfdataset(filesin,concat_dim='time')

    fin.coords['y']=fin.lat.values[0,:,0]
    fin.coords['x']=fin.lon.values[0,0,:]
    if var == 'PR':
        fin.PR.values = fin.PR.values*3600.
        fin.PR.attrs['units']='mm -hr-1'

    return fin

def load_obs_pr(obs_data,freq):

    path_obs = oinfo.get_param_values(obs_data,'path_obs')
    patt_obs = oinfo.get_param_values(obs_data,'patt_obs')
    filesin_obs = sorted(glob('%s/%s%s*' %(path_obs,patt_obs,freq)))
    print (filesin_obs)
    finobs = xr.open_mfdataset(filesin_obs)

    latvar = oinfo.get_param_values(obs_data,'latvar')
    lonvar = oinfo.get_param_values(obs_data,'lonvar')
    pr_varname = oinfo.get_param_values(obs_data,'pr_varname')
    scale = oinfo.get_param_values(obs_data,'scale')
    pr_varname = oinfo.get_param_values(obs_data,'pr_varname')
    scale = oinfo.get_param_values(obs_data,'scale')

    finobs.rename({latvar:'lat',lonvar:'lon'},inplace=True)
    finobs.rename({pr_varname:'PR'},inplace=True)


    finobs.PR.values = finobs.PR.values*scale
    finobs.PR.attrs['units']='mm hr-1'



    return finobs

def load_sel_seasons_wrf(wrun,freq,var,levs=[]):

    fin = load_wrf(wrun,freq,var,levs)
    finwrf_years = fin.sel(time=slice("%s-%s" %(cfg.syear,cfg.smonth),"%s-%s" %(cfg.eyear,cfg.emonth)))
    finwrf_seas  = finwrf_years.sel(time=is_month(finwrf_years['time.month']))

    return finwrf_seas

def load_sel_seasons_obs(obs,freq):

    finobs = load_obs_pr(obs,freq)
    finobs_years = finobs.sel(time=slice("%s-%s" %(cfg.syear,cfg.smonth),"%s-%s" %(cfg.eyear,cfg.emonth)))
    finobs_seas  = finobs_years.sel(time=is_month(finobs_years['time.month']))

    return finobs_seas

def load_sel_seasons_era5(fin):

    fin_years = fin.sel(time=slice("%s-%s" %(cfg.syear,cfg.smonth),"%s-%s" %(cfg.eyear,cfg.emonth)))
    fin_seas  = fin_years.sel(time=is_month(fin_years['time.month']))
    fin_seas.rename({'tp':'PR'},inplace=True)

    fin_seas.PR.values = fin_seas.PR.values*1000/24.

    return fin_seas

def load_sounding_YMC(path_sounding,dates_sounding,plevs):

    obs_all=pd.DataFrame()
    spl = np.zeros(len(plevs),dtype=np.float)

    for n, nom_date_sounding in enumerate(dates_sounding):
        date_filename = nom_date_sounding.strftime('%y%m%d%H')
        filein_sounding = glob("%s/??_%s.asc" %(path_sounding,date_filename))
        if len(filein_sounding)>1:
            print (filein_sounding)
            raise ValueError("Too many available files meeting date condition")
        obs = pd.read_csv(filein_sounding[0],delim_whitespace=True,header=11,skiprows=[12,13])

        for pl in range(len(plevs)):
            spl[pl]=obs.iloc[(obs['Pres']-plevs[pl]).abs().argsort()[:2]].Pres.values[0]

        obs_sel=obs[obs.Pres.isin(spl)].reset_index()
        obs_all=pd.concat((obs_all,obs_sel))

    obs_mean=obs_all.groupby(level=0).mean()

    slat = obs_mean.Lat.values[0]
    slon = obs_mean.Lon.values[0]

    #Calculating wind direction and speed (knts) from sounding

    sknt = ((obs_mean.Ucmp.values)**2+(obs_mean.Vcmp.values)**2)**(0.5)*1.94384
    drct = np.arctan2(obs_mean.Vcmp.values,obs_mean.Ucmp.values)*180./np.pi
    drct[drct<0]=360+drct[drct<0]

    sounding=dict(zip(('slat','slon','hght','pres','temp','dwpt','relh','drct','sknt'),(slat,slon,obs_mean.Alt.values,obs_mean.Pres.values,obs_mean.Temp.values,obs_mean.Dewpt.values,obs_mean.RH.values,drct,sknt)))

    return sounding

def load_sounding_WRF(wrun,dates_sounding,plevs,slat,slon):

    model_all=pd.DataFrame()
    fileref = nc.Dataset('%s/%s/out/wrfout_d01_2015-11-01_00:00:00' %(cfg.path_wrfout,wrun))
    varnames = ['TC','TD','RH','UA','VA','Z']

    for vname in varnames:
        for n, nom_date_sounding in enumerate(dates_sounding):
            date_filename = nom_date_sounding.strftime('%y%m%d%H')
            print('%s/%s/2015-2016/UIB_PLEVS_01H_%s_%s-%s.nc' %(cfg.path_in,wrun,vname,nom_date_sounding.year,nom_date_sounding.month))
            filesin = sorted(glob("%s/%s/2015-2016/UIB_PLEVS_01H_%s_%s-%s.nc" %(cfg.path_in,wrun,vname,nom_date_sounding.year,nom_date_sounding.month)))
            fin = xr.open_mfdataset(filesin)


            if len(filesin)==1:
                fin.coords['y']=fin.lat.values[:,0]
                fin.coords['x']=fin.lon.values[0,:]
            else:
                fin.coords['y']=fin.lat.values[0,:,0]
                fin.coords['x']=fin.lon.values[0,0,:]

            model= fin.sel(method='nearest',time=nom_date_sounding.strftime('%Y-%m-%d %H'),y=slat,x=slon)[vname].to_dataframe()
            model_all=pd.concat((model_all,model))

        model_mean=model_all.groupby(level=0).mean()


    model_mean['sknt'] = ((model_mean['UA'])**2+(model_mean['VA'])**2)**(0.5)*1.94384
    model_mean['drct'] = np.arctan2(model_mean['VA'],model_mean['UA'])*180./np.pi

    sounding=dict(zip(('hght','pres','temp','dwpt','relh','drct','sknt'),(model_mean['Z'],fin.levels.values,model_mean['TC'],model_mean['TD'],model_mean['RH'],model_mean['drct'],model_mean['sknt'])))

    return sounding
