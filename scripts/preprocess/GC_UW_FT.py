import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d
import math

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import ScalarMappable
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import function as f

# --- User inputs ---
year = 2019
month = 6

path = '../../data/GCdata/GC_output/'
file = '/GEOSChem.SpeciesConcHF.'+str(year)+str(month).zfill(2)+'01_0000z.nc4'
datapath = '../../data/GCdata/GC_processed/'

export = True
flux_export = True

bronx = [-73.8782,40.8680]
nbrunswick = [-74.4226,40.4728]
underhill = [-72.8684,44.5283]
huntington = [-74.2231,43.9731]

# --- end user inputs ---


'''
Define important functions:

conc_partition: Partitions FT and PBL based on GC-estimated PBLH
conc_partition_nest: do the same thing as conc_partition, except on a nested grid
upwind_concentration: defines the region that will be used to calculate upwind concentration from lat/lon bounds and level of pblh
species_sum: sums across all species (for instance, all wet deposition components)
NetFlux_Hg0: sums across all Hg0 flux inventories for a given lat/lon. This includes emissions and deposition!
'''

def conc_partition(conc,pblp,convert=True,upper_thresh=700,export=True,path='/PROVIDEPATH/',name='rural'):
    '''
    Script partitions free tropospheric and PBL concentrations based on GC values.
    Before running, make sure that conc has been read in fully (i.e., conc = conc.load())
    '''
    for i in list(range(0,len(conc.time))): 
        temp_r = conc.isel(time=i)
        thresh = pblp[i]

        pbl_all = temp_r['SpeciesConcVV_Hg0'].where(temp_r['P_m']>thresh,drop=True)
        check = math.isnan(pbl_all.mean().values.tolist())
        if check == True:
            pbl_i = temp_r['SpeciesConcVV_Hg0'].isel(lev=0)
        if check == False:
            pbl_i = pbl_all.mean('lev')
        ft_all = temp_r['SpeciesConcVV_Hg0'].where(temp_r['P_m']<thresh,drop=True).where(temp_r['P_m']>upper_thresh,drop=True)
        check = math.isnan(ft_all.mean().values.tolist())
        if check == True:
            ft_i = temp_r['SpeciesConcVV_Hg0'].sel(lev=upper_thresh,method='nearest')
        if check == False:
            ft_i = ft_all.mean('lev')

        if i == 0:
            pbl = pbl_i
            ft = ft_i
            levs = np.array(len(pbl_all))
            levs_ft = np.array(len(ft_all))
        else:
            pbl = xr.concat([pbl,pbl_i],dim='time')
            ft = xr.concat([ft,ft_i],dim='time')
            levs = np.append(levs,len(pbl_all))
            levs_ft = np.append(levs_ft,len(ft_all))
            
    if convert==True:       
        pbl = f.HgConversion(pbl)
        ft = f.HgConversion(ft)
    
    if export==True:
        np.save(path+name+'_pbl.npy',pbl)
        np.save(path+name+'_ft.npy',ft)
        np.save(path+name+'_levs_pbl.npy',levs)
        np.save(path+name+'_levs_ft.npy',levs_ft)

    return pbl,ft,levs,levs_ft


def conc_partition_nest(gc,pbl,ft,path='/PROVIDEPATH/',name='rural_nest',export=True):
    for i in list(range(0,len(pbl))):
        pbl_i = f.HgConversion(gc.isel(time=i,lev=slice(0,pbl[i]))['SpeciesConcVV_Hg0'].mean().values)
        ft_i = f.HgConversion(gc.isel(time=i,lev=slice(pbl[i],ft[i]))['SpeciesConcVV_Hg0'].mean().values)
        
        if i == 0:
            pbl_conc = np.array(pbl_i)
            ft_conc = np.array(ft_i)
        else:
            pbl_conc = np.append(pbl_conc,pbl_i)
            ft_conc = np.append(ft_conc,ft_i)
    if export == True:
        np.save(path+name+'_pbl.npy',pbl_conc)
        np.save(path+name+'_ft.npy',ft_conc)
        
    return pbl,ft

def upwind_concentration(gc,pbl,lonmin=-80,lonmax=-70,latmin=40,latmax=50,path='/PROVIDEPATH/',name='rural',export=True):
    for i in list(range(0,len(pbl))):
        uw_i = f.HgConversion(gc.isel(time=i,lev=slice(0,pbl[i])).sel(lat=slice(latmin,latmax),lon=slice(lonmin,lonmax))['SpeciesConcVV_Hg0'].mean().values)
        
        if i == 0:
            uw = np.array(uw_i)
        else:
            uw = np.append(uw,uw_i)
    
    if export == True:
        np.save(path+name+'_uw.npy',uw)
    return uw

def species_sum(ds,autolist=True,varlist=['WetLossLS_HgCl2']):
    #each element of ds must have the same units
    if autolist==True:
        for i in list(ds):
            if i == list(ds)[0]:
                ds_tot = ds[i]
            else:
                ds_tot = ds_tot+ds[i]
    if autolist==False:
        for i in varlist:
            if i == varlist[0]:
                ds_tot = ds[i]
            else:
                ds_tot = ds_tot+ds[i]
    return ds_tot

def NetFlux_Hg0(path,run,y='2014',m='06',savefig=True,hg0_only=False):

    drydep = xr.open_mfdataset(path+run+'/GEOSChem.DryDep.'+y+m+'01_0000z.nc4')
    ems = xr.open_mfdataset(path+run+'/GEOSChem.MercuryEmis.'+y+m+'01_0000z.nc4')
    hemco = xr.open_mfdataset(path+run+'/HEMCO_diagnostics.'+y+m+'01*0000.nc')

    drydep['DryDep_Hg0_ngm2hr'] = drydep['DryDep_Hg0']*100*100*(1/6.022e23)*200.59*1e9*3600
    
    ems['EmisHg0soil_ngm2hr'] = ems['EmisHg0soil']/drydep['AREA']*1e12*3600
    ems['EmisHg0ocean_ngm2hr'] = ems['EmisHg0ocean']/drydep['AREA']*1e12*3600
    ems['EmisHg0land_ngm2hr'] = ems['EmisHg0land']/drydep['AREA']*1e12*3600
    ems['EmisHg0snow_ngm2hr'] = ems['EmisHg0snow']/drydep['AREA']*1e12*3600
    
    hemco['EmisHg0_Natural_ngm2hr'] = hemco['EmisHg0_Natural']*1e12*3600
    hemco['EmisHg0_BioBurn_ngm2hr'] = hemco['EmisHg0_BioBurn']*1e12*3600
    hemco['EmisHg0_ASGM_ngm2hr'] = hemco['EmisHg0_ASGM']*1e12*3600
    hemco['EmisHg0_Anthro_ngm2hr'] = hemco['EmisHg0_Anthro']*1e12*3600
    
    drydep_avg  = drydep.groupby('time.month').mean()[['AREA','DryDep_Hg0_ngm2hr']]
    ems_avg     = ems.groupby('time.month').mean()[['EmisHg0snow_ngm2hr','EmisHg0land_ngm2hr',"EmisHg0soil_ngm2hr","EmisHg0ocean_ngm2hr"]]
    hemco_avg   = hemco.groupby('time.month').mean()[['EmisHg0_ASGM_ngm2hr','EmisHg0_Anthro_ngm2hr',"EmisHg0_BioBurn_ngm2hr","EmisHg0_Natural_ngm2hr"]]

    averagefluxes = xr.merge([drydep_avg,
                              ems_avg,
                              hemco_avg])
    return averagefluxes

def geos_chem_n_parameters(path,file,datapath,rural_levs,rural_levs_ft,urban_levs,urban_levs_ft,
                           run='run0051',
                           export=True,
                           method='midpoint',
                           lons_r=[-76.6,-74.65],
                           lats_r=[42.25,43.25],
                           lons_u=[-76.6,-74.65],
                           lats_u=[39.75,40.75]):
    
    ds = xr.open_mfdataset(path+run+file).load()
    
    if method == 'midpoint':
        r_lon = (huntington[0]+underhill[0])/2
        r_lat = (huntington[1]+underhill[1])/2
        u_lon = (nbrunswick[0]+bronx[0])/2
        u_lat = (nbrunswick[1]+bronx[1])/2
    else:
        r_lon = huntington[0]
        r_lat = huntington[1]
        u_lon = nbrunswick[0]
        u_lat = nbrunswick[1]
        
    rural = ds.sel(lon=r_lon,lat=r_lat,method='nearest')
    urban = ds.sel(lon=u_lon,lat=u_lat,method='nearest')
            
    #partition concentrations at the selected gridcell into pbl and ft, export pbl and ft concentrations
    rural_conc_pbl_n,rural_conc_ft_n = conc_partition_nest(rural,
                                                rural_levs,
                                                rural_levs_ft,
                                                export=export,
                                                path=datapath+str(month).zfill(2)+'/',
                                                name='rural_n_'+run+'_'+str(round(r_lat,2))+'_'+str(round(r_lon,2)))
    urban_conc_pbl_n,urban_conc_ft_n = conc_partition_nest(urban,
                                                urban_levs,
                                                urban_levs_ft,
                                                export=export,
                                                path=datapath+str(month).zfill(2)+'/',
                                                name='urban_n_'+run+'_'+str(round(u_lat,2))+'_'+str(round(u_lon,2)))
    
    #calculating uw concentrations
    rural_uw_n = upwind_concentration(ds,
                                      rural_levs,
                                      lonmin=lons_r[0],
                                      lonmax=lons_r[1],
                                      latmin=lats_r[0],
                                      latmax=lats_r[1],
                                      path=datapath+str(month).zfill(2)+'/',
                                      name='rural_n_'+run+'_'+str(lats_r[0])+'_'+str(lats_r[1])+'_'+str(lons_r[0])+'_'+str(lons_r[1]),
                                      export=export)
    urban_uw_n = upwind_concentration(ds,
                                      urban_levs,
                                      lonmin=lons_u[0],
                                      lonmax=lons_u[1],
                                      latmin=lats_u[0],
                                      latmax=lats_u[1],
                                      path=datapath+str(month).zfill(2)+'/',
                                      name='urban_n_'+run+'_'+str(lats_u[0])+'_'+str(lats_u[1])+'_'+str(lons_u[0])+'_'+str(lons_u[1]),
                                      export=export)
    if export == False:
        return rural_conc_pbl_n,rural_conc_ft_n,rural_uw_n,urban_conc_pbl_n,urban_conc_ft_n,urban_uw_n


#writing single function that will do all of the things above for a given file, export files using desired names.
#output file names should specify the run they are based on, and the horizontal and vertical thresholds used.
#not wrapping the reading of standard res run, info needed for nested calculations.

#defining the 1976 std atmosphere (using values from S&P)
p_stdatm = np.array([1013,898,795,701.2,616.6,540.5,472.2,411.1,356.5,308,265,227,194,165,142,121,103,88.5,75.6,64.6,5.53]) #converts from hPa to Pa
T_stdatm = np.array([288,282,275,269,262,256,249,243,236,230,223,217,217,217,217,217,217,217,217,217,217]) #K
z_stdatm = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])*1000 #m

p_to_z = interp1d(p_stdatm,z_stdatm,fill_value='extrapolate')
z_to_p = interp1d(z_stdatm,p_stdatm,fill_value='extrapolate')

#reading in necessary met (using 2x2.5 for efficiency, nested makes little difference when aggregated across region)
met = xr.open_dataset('../../data/pblh/MERRA2/MERRA2.'+str(year)+str(month).zfill(2)+'.pblh.2x25.nc4')

#now reading in concentrations (edit this based on your needs
base_s = xr.open_mfdataset(path+'STND'+file).load()

#calculating pressure at each gridbox using standard atmosphere and definition of GEOS-Chem hybrid-sigma coords
p_i = base_s['hyai']+(base_s['hybi']*(met['SLP'])/100)
p_m = p_i.rolling(ilev=2).mean().dropna('ilev')
p_m = p_m.assign_coords(ilev=(base_s.lev.values))
p_m = p_m.rename({'ilev':'lev'})
p_m = p_m.sel(lat=base_s['lat'],lon=base_s['lon']).resample(time='1H').mean()
base_s = xr.merge([base_s,p_m.to_dataset(name='P_m')])

pblh_rural = met['PBLH'].sel(lon=huntington[0],lat=huntington[1],method='nearest')
pblh_urban = met['PBLH'].sel(lon=nbrunswick[0],lat=nbrunswick[1],method='nearest')

pblp_rural = z_to_p(pblh_rural.values)
pblp_urban = z_to_p(pblh_urban.values)

rural = base_s.sel(lon=huntington[0],lat=huntington[1],method='nearest')
urban = base_s.sel(lon=nbrunswick[0],lat=nbrunswick[1],method='nearest')


#getting boundary concentrations for 2x2.5 run
rural_conc_pbl,rural_conc_ft,rural_levs,rural_levs_ft = conc_partition(rural,
                                                                       pblp_rural,
                                                                       export=export,
                                                                       path=datapath,
                                                                       name='rural_'+str(month).zfill(2))
urban_conc_pbl,urban_conc_ft,urban_levs,urban_levs_ft = conc_partition(urban,
                                                                       pblp_urban,
                                                                       export=export,
                                                                       path=datapath,
                                                                       name='urban_'+str(month).zfill(2))

rural_uw_s = upwind_concentration(base_s,
                                  rural_levs,
                                  lonmin=-76.6,
                                  lonmax=-74.65,
                                  latmin=42.25,
                                  latmax=44.25,
                                  path=datapath,
                                  name='rural_'+str(month).zfill(2))
urban_uw_s = upwind_concentration(base_s,
                                  urban_levs,
                                  lonmin=-76.6,
                                  lonmax=-74.65,
                                  latmin=39.75,
                                  latmax=40.75,
                                  path=datapath,
                                  name='urban_'+str(month).zfill(2))

# --- exporting nested variables ---
# --- std ---
print('starting std')
geos_chem_n_parameters(path,
                     file,
                     datapath,
                     rural_levs,
                     rural_levs_ft,
                     urban_levs,
                     urban_levs_ft,
                     run='NEST',
                     export=export,
                     method='midpoint')

# --- 3x ---
print('starting 3x')
geos_chem_n_parameters(path,
                     file,
                     datapath,
                     rural_levs,
                     rural_levs_ft,
                     urban_levs,
                     urban_levs_ft,
                     run='ANT3-ENA',
                     export=export,
                     method='midpoint')

# --- 10x ---
print('starting 10x')
geos_chem_n_parameters(path,
                     file,
                     datapath,
                     rural_levs,
                     rural_levs_ft,
                     urban_levs,
                     urban_levs_ft,
                     run='ANT10-ENA',
                     export=export,
                     method='midpoint')

# --- 7x ---
print('starting 7x')
geos_chem_n_parameters(path,
                     file,
                     datapath,
                     rural_levs,
                     rural_levs_ft,
                     urban_levs,
                     urban_levs_ft,
                     run='ANT7-ENA',
                     export=export,
                     method='midpoint')
    
# --- 3x_NEscale ---
print('starting 3x_NEscale')
geos_chem_n_parameters(path,
                     file,
                     datapath,
                     rural_levs,
                     rural_levs_ft,
                     urban_levs,
                     urban_levs_ft,
                     run='ANT3',
                     export=export,
                     method='midpoint')

# --- 10x_NEscale ---
print('starting 10x_NEscale')
geos_chem_n_parameters(path,
                     file,
                     datapath,
                     rural_levs,
                     rural_levs_ft,
                     urban_levs,
                     urban_levs_ft,
                     run='ANT10',
                     export=export,
                     method='midpoint')

# --- 7x_NEscale ---
print('starting 7x_NEscale')
geos_chem_n_parameters(path,
                     file,
                     datapath,
                     rural_levs,
                     rural_levs_ft,
                     urban_levs,
                     urban_levs_ft,
                     run='ANT7',
                     export=export,
                     method='midpoint')
    

# --- Free tropospheric scaling sensitivity ---
print('starting scaled FT concentrations')
geos_chem_n_parameters(path,
                     file,
                     datapath,
                     rural_levs,
                     rural_levs_ft,
                     urban_levs,
                     urban_levs_ft,
                     run='CFTSCALMAX',
                     export=export,
                     method='midpoint')

# --- Free tropospheric scaling sensitivity (tuned) ---
print('starting scaled FT concentrations (tuned)')
geos_chem_n_parameters(path,
                     file,
                     datapath,
                     rural_levs,
                     rural_levs_ft,
                     urban_levs,
                     urban_levs_ft,
                     run='CFTSCAL',
                     export=export,
                     method='midpoint')


if flux_export == True:
    print('calculating net fluxes')
    emslist = ['EmisHg0snow_ngm2hr',
               'EmisHg0land_ngm2hr',
               'EmisHg0soil_ngm2hr',
               'EmisHg0ocean_ngm2hr',
               'EmisHg0_ASGM_ngm2hr',
               'EmisHg0_Anthro_ngm2hr',
               'EmisHg0_BioBurn_ngm2hr',
               'EmisHg0_Natural_ngm2hr']

    # --- std ---
    run = 'NEST'
    name = 'NetHg0Flux_'

    averagefluxes = NetFlux_Hg0(path,run,y=str(year),m=str(month).zfill(2))
    averagefluxes['EmisTot_ngm2hr'] = species_sum(averagefluxes,autolist=False,varlist=emslist)
    averagefluxes.to_netcdf(path=datapath+name+str(month).zfill(2)+'.nc4')

    # --- 3x ---
    run = 'ANT3-ENA'
    name = 'NetHg0Flux_3x_'

    averagefluxes = NetFlux_Hg0(path,run,y=str(year),m=str(month).zfill(2))
    averagefluxes['EmisTot_ngm2hr'] = species_sum(averagefluxes,autolist=False,varlist=emslist)
    averagefluxes.to_netcdf(path=datapath+name+str(month).zfill(2)+'.nc4')

    # --- 10x ---
    run = 'ANT10-ENA'
    name = 'NetHg0Flux_10x_'

    averagefluxes = NetFlux_Hg0(path,run,y=str(year),m=str(month).zfill(2))
    averagefluxes['EmisTot_ngm2hr'] = species_sum(averagefluxes,autolist=False,varlist=emslist)
    averagefluxes.to_netcdf(path=datapath+name+str(month).zfill(2)+'.nc4')
    
    # --- 7x ---
    run = 'ANT7-ENA'
    name = 'NetHg0Flux_7x_'

    averagefluxes = NetFlux_Hg0(path,run,y=str(year),m=str(month).zfill(2))
    averagefluxes['EmisTot_ngm2hr'] = species_sum(averagefluxes,autolist=False,varlist=emslist)
    averagefluxes.to_netcdf(path=datapath+name+str(month).zfill(2)+'.nc4')
    
    # --- 3x_NEscale ---
    run = 'ANT3'
    name = 'NetHg0Flux_3x_NEscale_'

    averagefluxes = NetFlux_Hg0(path,run,y=str(year),m=str(month).zfill(2))
    averagefluxes['EmisTot_ngm2hr'] = species_sum(averagefluxes,autolist=False,varlist=emslist)
    averagefluxes.to_netcdf(path=datapath+name+str(month).zfill(2)+'.nc4')

    # --- 10x ---
    run = 'ANT10'
    name = 'NetHg0Flux_10x_NEscale_'

    averagefluxes = NetFlux_Hg0(path,run,y=str(year),m=str(month).zfill(2))
    averagefluxes['EmisTot_ngm2hr'] = species_sum(averagefluxes,autolist=False,varlist=emslist)
    averagefluxes.to_netcdf(path=datapath+name+str(month).zfill(2)+'.nc4')
    
    # --- 7x ---
    run = 'ANT7'
    name = 'NetHg0Flux_7x_NEscale_'

    averagefluxes = NetFlux_Hg0(path,run,y=str(year),m=str(month).zfill(2))
    averagefluxes['EmisTot_ngm2hr'] = species_sum(averagefluxes,autolist=False,varlist=emslist)
    averagefluxes.to_netcdf(path=datapath+name+str(month).zfill(2)+'.nc4')

    # --- scaled FT concentrations ---
    # --- NOTE: total emissions need factor of 1e-14 kg m-2 s-1 subtracted - this is the magnitude of emissions added at 3500 m to scale FT conc. ---
    run = 'CFTSCALMAX'
    name = 'NetHg0Flux_FTscale_'

    averagefluxes = NetFlux_Hg0(path,run,y=str(year),m=str(month).zfill(2))
    averagefluxes['EmisTot_ngm2hr'] = species_sum(averagefluxes,autolist=False,varlist=emslist)
    averagefluxes.to_netcdf(path=datapath+name+str(month).zfill(2)+'.nc4')

    # --- scaled FT concentrations (tuned) ---
    # --- NOTE: total emissions need factor of 5e-15 kg m-2 s-1 subtracted - this is the magnitude of emissions added at 3500 m to scale FT conc. ---
    run = 'CFTSCAL'
    name = 'NetHg0Flux_FTscale_tuned_'

    averagefluxes = NetFlux_Hg0(path,run,y=str(year),m=str(month).zfill(2))
    averagefluxes['EmisTot_ngm2hr'] = species_sum(averagefluxes,autolist=False,varlist=emslist)
    averagefluxes.to_netcdf(path=datapath+name+str(month).zfill(2)+'.nc4')


    

