

#Import modules
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d
from scipy.stats import norm
from scipy.stats import uniform
import function as f

#user inputs
year = 2017
month = 12
option = 'mean'  #optimization strategy
roll = 5
wthresh = 25

start = '2014-01-01' #start of considered obs
end = '2020-01-01'  #end of considered obs

uw = np.linspace(1,2,40)
ft = np.linspace(1,2,40)
fluxes = np.linspace(-75,75,200)

amnetpath = '/net/fs03/d0/emroy/paper_staging/RegionalBoxModel/data/NADP/AMNET-NE-h.csv'
outpath   = '/net/fs03/d0/emroy/paper_staging/RegionalBoxModel/data/boxmod_opt/'

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


flux_urban_array = np.array([])
flux_rural_array = np.array([])

for i in ft:
    fua = np.array([])
    fra = np.array([])
    for j in uw:

        flux_urban,junk = f.optimal_flux(pbl,j,i,fluxes,measurements,var='urban',option=option)
        flux_rural,junk = f.optimal_flux(pbl,j,i,fluxes,measurements,var='rural',option=option)
        
        fua = np.append(fua,flux_urban)
        fra = np.append(fra,flux_rural)

    if i == ft[0]:
        flux_urban_array = fua
        flux_rural_array = fra
    else:
        flux_urban_array = np.vstack((flux_urban_array,fua))
        flux_rural_array = np.vstack((flux_rural_array,fra))

#now exporting both 2d np arrays
np.save(outpath+'flux_urban_solspace_'+str(month).zfill(2)+'.npy',flux_urban_array)
np.save(outpath+'flux_rural_solspace_'+str(month).zfill(2)+'.npy',flux_rural_array)
