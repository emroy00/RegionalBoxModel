"""
June 30, 2025

All box modeling functions will live in this script. Here, they will be callable to separate functions that evaluate across predefined solution spaces and produce key figures for the paper.
"""

#Import modules
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d
from scipy.stats import norm
from scipy.stats import uniform


#defining important functions

def HgConversion(conc):
    '''
    Converts Hg concentration to ng m3
    '''
    
    mw_hg = 200.59 #g/mol
    std_p = 101325 #Pa
    std_T = 273.15 #K
    R = 8.314 #m3 Pa K-1 mol-1
    ng = 1e9
    ngm3 = conc*mw_hg*std_p*ng/(R*std_T)
    return ngm3

def pblh_read(path,tz_convert=False,tz='est'):
    '''
    Opens relevant pblh dataset from Zhang et al. (2019) and converts to a pandas dataframe.
    '''
    ds = xr.open_dataset(path)
    
    ds['year'] = ds.time[0]
    ds['month'] = ds.time[1]
    ds['day'] = ds.time[2]
    ds['hour'] = ds.time[3]
    ds['datetime64'] = pd.to_datetime(pd.DataFrame({'year':ds['year'].values,'month':ds['month'].values,'day':ds['day'].values,'hour':ds['hour']}))

    ds_out = xr.Dataset(
        data_vars=dict(
            pblh=(["time"], ds['pblh'].values),
            u=(["height","time"], ds['u'].values),
            v=(["height","time"], ds['v'].values),
            T=(["height","time"], ds['t'].values)),
        coords=dict(
            time=ds['datetime64'].values,
            height=ds['h_magl'].values))
    if tz_convert==True:
        ds_out = ds_out.to_dataframe()#.tz_localize(tz='UTC').tz_convert('America/New_York').tz_localize(None)

    return ds_out


def open_pbl(loc='EWR',year=2017,month=6,lst='EST',roll=5,wthresh=100,dirpath='../../data/pblh/'):
    pbl = pblh_read(dirpath+loc+'_20m_interp_profiles_xr.nc')
    pbl = pbl.sel(time=str(year)+'-'+str(month).zfill(2))

    #if only considering EST (NADP reporting policy)
    if lst == 'EST':
        shift = -5
    if lst == 'EDT':
        shift = -4

    #if switching between EDT and EST
    #if month in range(4, 11):
    #    shift=4
    #else:
    #    shift=5

    pbl = pbl.shift(time=shift) #shifting from UTC to EST, note that sign of change is different from the commented approach below
    
    #pbl['time'] = pbl.time.shift(time=shift) #shifting from UTC to EST, note that sign of change is different from accepted approach above.
    
    pbl = pbl.sel(time=str(year)+'-'+str(month).zfill(2))
    pbl['u'] = pbl['u'].where(np.abs(pbl['u'])<wthresh)
    pbl['v'] = pbl['v'].where(np.abs(pbl['v'])<wthresh)
    pbl['hour'] = pbl.time.dt.hour
    pbl = pbl.groupby('hour').mean()
    
    windmag = np.array([])
    for i in pbl.hour.values:
        pbl_height = pbl.sel(hour=i)['pblh'].values
        in_pbl = pbl.sel(hour=i,height=slice(0,pbl_height))
        in_pbl = (in_pbl['u']**2+in_pbl['v']**2)**0.5
        in_pbl = in_pbl.mean()
        windmag = np.append(windmag,in_pbl.values)

    pbl_df = pbl['pblh'].to_dataframe()
    pbl_df['pblh_roll'] = pbl_df['pblh'].rolling(roll,center=True,min_periods=1).mean()
    pbl_df['windmag'] = windmag
    pbl_df['windmag_roll'] = pbl_df['windmag'].rolling(roll,center=True,min_periods=1).mean()
    
    return pbl, pbl_df


def box_model(flux,pbl,measurements,dt=3600,dy=100000,start_date='6/1/2017',end_date='6/15/2017',spinup=5,groupby=False,obsmerge=True,upw_conc=1,hg_ft=1.5,pbl_var='avg_pblh_roll',wind_var='avg_windmag_roll',roll=5):
    
    '''
    Based on Seinfeld and Pandis Ch 25, equations 25.13 and 25.14
    '''
    
    time = pd.date_range(start=start_date,end=end_date,freq='1H')
    z_s = 0                                                                            #m
    m_air = 0.02897                                                                    #kg mol-1

    hg_box = 2                 #ng m-3

    box_ar = np.array(hg_box)
    ent_ar = np.array([])
    upw_ar = np.array([])
    flx_ar = np.array([])
    
    for i in list(range(0,len(time))):
        pbl_index = int(np.rint((i/24-np.floor(i/24))*24)) #This provides the number of the hour for each day of the timeseries, will take from average pbl
        z_pbl = pbl.iloc[pbl_index][pbl_var]               #m
        z_pblm1 = pbl.iloc[pbl_index-1][pbl_var]           #m
        windmag = pbl.iloc[pbl_index][wind_var]            #m s-1

        dpbl_z = z_pbl-z_pblm1                             #z hr-1

        if dpbl_z>0:
            hg_ent = ((dpbl_z/(z_pbl))*(hg_ft-hg_box))     #ng m-3 hr-1
        else:
            hg_ent = 0

        hg_upw = (windmag/dy)*(upw_conc-hg_box)*dt         #ng m-3 hr-1; windmag in units of m s-1, dt=3600 converts to m hr-1

        hg_flx = flux/z_pbl                                #ng m-3 hr-1; flux in units of ng m-2 hr-1

        hg_box = hg_box+hg_ent+hg_upw+hg_flx               #hg_box represents new concentration after 1 hour timestep

        box_ar = np.append(box_ar,hg_box)
        ent_ar = np.append(ent_ar,hg_ent)
        upw_ar = np.append(upw_ar,hg_upw)
        flx_ar = np.append(flx_ar,hg_flx)
    boxmodel = pd.DataFrame({'time':time,'hg_box':box_ar[1:],'hg_ent':ent_ar,'hg_upw':upw_ar,'hg_flx':flx_ar})
    boxmodel = boxmodel.set_index('time')
    boxmodel = boxmodel[boxmodel.index>time[24*spinup]]
    boxmodel['hour'] = boxmodel.index.hour
    
    if groupby == True:
        boxmodel = boxmodel.groupby('hour').mean()
        if obsmerge == True:
            boxmodel = pd.merge(boxmodel,measurements,left_index=True,right_index=True)
            
            #boxmodel = pd.merge(boxmodel,NY20['GEM'],left_index=True,right_index=True)
            #boxmodel = boxmodel.rename(columns={'GEM':'NY20'})
            #boxmodel = pd.merge(boxmodel,NY06['GEM'],left_index=True,right_index=True)
            #boxmodel = boxmodel.rename(columns={'GEM':'NY06'})
            #boxmodel = pd.merge(boxmodel,NJ30['GEM'],left_index=True,right_index=True)
            #boxmodel = boxmodel.rename(columns={'GEM':'NJ30'})
            #boxmodel = pd.merge(boxmodel,VT99['GEM'],left_index=True,right_index=True)
            #boxmodel = boxmodel.rename(columns={'GEM':'VT99'})
            #boxmodel['urban'] = ((boxmodel['NY06'].interpolate()+boxmodel['NJ30'].interpolate())/2).rolling(roll,center=True,min_periods=1).mean()
            #boxmodel['rural'] = ((boxmodel['NY20'].interpolate()+boxmodel['VT99'].interpolate())/2).rolling(roll,center=True,min_periods=1).mean()
        return boxmodel
    else:
        return boxmodel

    
def optimal_flux(obs,upw_conc,hg_ft,fluxes,measurements,var='urban',option='standard'):

    prob = 0
    f_opt = 0
    measure = measurements[var]
    #measure = box_model(0,obs,measurements,upw_conc=upw_conc,hg_ft=hg_ft,groupby=True,obsmerge=True)[var]

    for f in fluxes:
        mod = box_model(f,obs,measurements,upw_conc=upw_conc,hg_ft=hg_ft,groupby=True,obsmerge=False)['hg_box']
        
        if option=='standard':
            mod_prob = norm.pdf(mod,loc=measure,scale=0.05).mean()
        if option=='mean':
            mod_prob = norm.pdf(mod.mean(),loc=measure.mean())

        if mod_prob>prob:
            f_opt = f
            prob = mod_prob

    return f_opt,prob

def zerofilter(data,thresh=1e-20):
    zf = data.where(data>thresh)    
    return zf

def amplitude_calc(series):
    return series.max()-series.min()

