#Import modules
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d
from scipy.stats import norm
from scipy.stats import uniform
from scipy.optimize import minimize
import function as f

#user inputs
year = 2017
month = 6
option = 'mean'  #optimization strategy
roll = 5
wthresh = 25

start = '2014-01-01' #start of considered obs
end = '2020-01-01'  #end of considered obs

uw_array = np.linspace(1,2,5)
ft_array = np.linspace(1,2,5)
fluxes = np.linspace(-75,75,200)

amnetpath = '/net/fs03/d0/emroy/paper_staging/RegionalBoxModel/data/NADP/AMNET-NE-h.csv'
outpath = '/net/fs03/d0/emroy/paper_staging/RegionalBoxModel/data/boxmod_opt/'

EWR,EWR_df = f.open_pbl(month=month,wthresh=wthresh)
BDL,BDL_df = f.open_pbl(month=month,loc='BDL',wthresh=wthresh)

### processing pbl data for nearby airports ###
pbl = pd.DataFrame({'EWR_pblh':EWR_df['pblh'],
                        'BDL_pblh':BDL_df['pblh'],
                        'EWR_windmag':EWR_df['windmag'],
                        'BDL_windmag':BDL_df['windmag']})

pbl['avg_pblh'] = pbl[['EWR_pblh', 'BDL_pblh']].mean(axis=1, skipna=True) #, 'MHT_pblh'
pbl['avg_pblh_roll'] = pbl['avg_pblh'].rolling(roll,center=True,min_periods=1).mean()

pbl['avg_windmag'] = pbl[['EWR_windmag', 'BDL_windmag']].mean(axis=1, skipna=True) #, 'MHT_windmag'
pbl['avg_windmag_roll'] = pbl['avg_windmag'].rolling(roll,center=True,min_periods=1).mean()


### Reading in NADP data ###

Amnet=pd.read_csv(amnetpath,na_values=(-9))
Amnet['collStart'] = pd.to_datetime(Amnet['collStart'])
Amnet['collEnd'] = pd.to_datetime(Amnet['collEnd'])
Amnet = Amnet.set_index('collStart')

#subsetting sites of interest
start = '2014-01-01'
end = '2020-01-01'

Amnet = Amnet[Amnet.index>start]
Amnet = Amnet[Amnet.index<end]
Amnet_all = Amnet.copy()
Amnet = Amnet[Amnet.index.month==month]
Amnet['hour'] = Amnet.index.hour 

NJ54 = Amnet[Amnet['SiteID']=='NJ54'].groupby('hour').mean()
NJ30 = Amnet[Amnet['SiteID']=='NJ30'].groupby('hour').mean()
NY06 = Amnet[Amnet['SiteID']=='NY06'].groupby('hour').mean()
NY20 = Amnet[Amnet['SiteID']=='NY20'].groupby('hour').mean()
VT99 = Amnet[Amnet['SiteID']=='VT99'].groupby('hour').mean()
MD08 = Amnet[Amnet['SiteID']=='MD08'].groupby('hour').mean()
MD99 = Amnet[Amnet['SiteID']=='MD98'].groupby('hour').mean()
NJ05 = Amnet[Amnet['SiteID']=='NJ05'].groupby('hour').mean()

NJ30_h = Amnet_all[Amnet_all['SiteID']=='NJ30']
NY06_h = Amnet_all[Amnet_all['SiteID']=='NY06']
NY20_h = Amnet_all[Amnet_all['SiteID']=='NY20']
VT99_h = Amnet_all[Amnet_all['SiteID']=='VT99']

#making dataframe containing all Hg0 observations and regional averages

measurements = pd.merge(NY20['GEM'],NY06['GEM'],left_index=True,right_index=True)
measurements = measurements.rename(columns={'GEM_x':'NY20','GEM_y':'NY06'})
measurements = pd.merge(measurements,NJ30['GEM'],left_index=True,right_index=True)
measurements = measurements.rename(columns={'GEM':'NJ30'})
measurements = pd.merge(measurements,VT99['GEM'],left_index=True,right_index=True)
measurements = measurements.rename(columns={'GEM':'VT99'})
measurements['urban'] = ((measurements['NY06'].interpolate()+measurements['NJ30'].interpolate())/2).rolling(roll,center=True,min_periods=1).mean()
measurements['rural'] = ((measurements['NY20'].interpolate()+measurements['VT99'].interpolate())/2).rolling(roll,center=True,min_periods=1).mean()

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

        if len(flux) == 1:
            hg_flx = flux/z_pbl                                #ng m-3 hr-1; flux in units of ng m-2 hr-1
        if len(flux) == 24:
            hg_flx = flux[pbl_index]/z_pbl

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
        return boxmodel
    else:
        return boxmodel

def bm_wrapper_urb(x):
    a,b,c,d,e,f = x
    
    p0 = np.ones(4)*a
    p1 = np.ones(4)*b
    p2 = np.ones(4)*c
    p3 = np.ones(4)*d
    p4 = np.ones(4)*e
    p5 = np.ones(4)*f
    
    flx_24 = np.append([p0],[p1,p2,p3,p4,p5])
    
    mod = box_model(flx_24,pbl,measurements,upw_conc=uw,hg_ft=ft,groupby=True,obsmerge=False)['hg_box']
    
    return sum((mod - measurements['urban'])**2)

def bm_wrapper_rur(x):
    a,b,c,d,e,f = x
    
    p0 = np.ones(4)*a
    p1 = np.ones(4)*b
    p2 = np.ones(4)*c
    p3 = np.ones(4)*d
    p4 = np.ones(4)*e
    p5 = np.ones(4)*f
    
    flx_24 = np.append([p0],[p1,p2,p3,p4,p5])
    
    mod = box_model(flx_24,pbl,measurements,upw_conc=uw,hg_ft=ft,groupby=True,obsmerge=False)['hg_box']
    
    return sum((mod - measurements['rural'])**2)


def bm_converter(x_opt):
    a,b,c,d,e,f = x_opt
    
    p0 = np.ones(4)*a
    p1 = np.ones(4)*b
    p2 = np.ones(4)*c
    p3 = np.ones(4)*d
    p4 = np.ones(4)*e
    p5 = np.ones(4)*f
    
    flx_24 = np.append([p0],[p1,p2,p3,p4,p5])
    
    mod = box_model(flx_24,pbl,measurements,upw_conc=uw,hg_ft=ft,groupby=True,obsmerge=False)['hg_box']
    
    return mod

for ft in ft_array:
    fua = np.array([])
    fra = np.array([])
    
    for uw in uw_array:

        #efficient approach
        flux_urban,junk = f.optimal_flux(pbl,uw,ft,fluxes,measurements,var='urban',option=option)
        flux_rural,junk = f.optimal_flux(pbl,uw,ft,fluxes,measurements,var='rural',option=option)
        
        fua = np.append(fua,flux_urban)
        fra = np.append(fra,flux_rural)
        
        #scipy approach
        x0 = np.ones(6)*0
        res_urb = minimize(bm_wrapper_urb, x0, tol=1e-3)
        
        x0 = np.ones(6)*0
        res_rur = minimize(bm_wrapper_rur, x0, tol=1e-3)

        if uw == uw_array[0]:
            fua_sci = res_urb.x
            fra_sci = res_rur.x
        else:
            fua_sci = np.vstack((fua_sci,res_urb.x))
            fra_sci = np.vstack((fra_sci,res_rur.x))
        
    if ft == ft_array[0]:
        flux_urban_array = fua
        flux_rural_array = fra
        
        flux_urban_array_sci = fua_sci
        flux_rural_array_sci = fra_sci
    else:
        flux_urban_array = np.vstack((flux_urban_array,fua))
        flux_rural_array = np.vstack((flux_rural_array,fra))
        
        flux_urban_array_sci = np.dstack((flux_urban_array_sci,fua_sci))
        flux_rural_array_sci = np.dstack((flux_rural_array_sci,fra_sci))

#now exporting both 2d np arrays
np.save(outpath+'flux_urban_solspace_supplemental_'+str(month).zfill(2)+'.npy',flux_urban_array)
np.save(outpath+'flux_rural_solspace_supplemental_'+str(month).zfill(2)+'.npy',flux_rural_array)

np.save(outpath+'flux_urban_solspace_supplementalscipy_'+str(month).zfill(2)+'.npy',flux_urban_array_sci) #NOTE: MUST BE TRANSPOSED BEFORE COMPARING MEANS
np.save(outpath+'flux_rural_solspace_supplementalscipy_'+str(month).zfill(2)+'.npy',flux_rural_array_sci) #NOTE: MUST BE TRANSPOSED BEFORE COMPARING MEANS
