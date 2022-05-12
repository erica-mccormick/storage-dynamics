import xarray
import scipy
import pandas as pd
import numpy as np
import os
import scipy.stats
from scipy import stats
import os.path
from os import path
import rioxarray
import rasterio as rio
import matplotlib.pyplot as plt
import matplotlib

# Specify scenario and models
scenario = 'rcp45'
model_names = ['CCSM4', 'CESM1-BGC','CMCC-CMS', 'CanESM2','GFDL-CM3','HadGEM2-CC', 'MIROC5','ACCESS1-0','HadGEM2-ES','CNRM-CM5']
models = [i + '_' + scenario for i in model_names]

def get_percentile_array():
    """
    Returns 10 arrays (one for each model) with the percentile rank of Smax
    relative to the future (2060-2100) annual precipitation distribution. 
    """
    # Import Smax and model data
    smax = rioxarray.open_rasterio('starting_data/Sr_2003_2020.tif').sel(band=1)
    dataset_path = 'starting_data/rcp45_models.nc'

    for model_run in models:
        p = xarray.open_dataset(dataset_path)[model_run].sel(year=slice(2060, 2100))
        p.coords['lon'] = p['lon'] - 360
        new_lon = smax.coords['x'].values
        new_lat = smax.coords['y'].values
        #pi = interpolated p
        pi = p.interp(coords = {"lat":new_lat, "lon":new_lon}, method = 'nearest')
        pi_mm = pi * 3.154e+7 # num of seconds in a year

        x_shape = len(pi_mm.coords['lon'].values)
        y_shape = len(pi_mm.coords['lat'].values)
        blank_array = np.empty([x_shape, y_shape])
        
        for x in range(x_shape):
            for y in range(y_shape):
                blank_array[x,y] = scipy.stats.percentileofscore(pi_mm.sel(lon = pi_mm.coords['lon'].values[x], lat = pi_mm.coords['lat'].values[y]).values, sr.sel(x = pi_mm.coords['lon'].values[x], y = pi_mm.coords['lat'].values[y]).values, kind = 'weak')

        with open('arrays/' + scenario + '_percentile_' + model_run + '.csv', 'w') as f:
            np.savetxt(f, blank_array)

        print('Done with model ' + model_run)
        
def write_arrays_to_tifs():
    """
    Convert arrays from get_percentile_array() to tifs. 
    """
    for model_runs in models:
        perc = np.loadtxt('arrays/'+scenario+'/'+ scenario + '_percentile_' + model_runs + '.csv')
        perc_2 = np.swapaxes(perc, 0, 1) # Necessary because of flipped x,y error in future_percentile.py

        # Get transformation data from Sr tif
        tif = rio.open('starting_data/Sr_2003_2020.tif')
        transform = tif.transform
        array_sr = tif.read(1)
        shape = array_sr.shape

        array_to_export = perc_2
        new_dataset = rio.open('tifs/'+scenario+'/percentile_' + model_runs +'.tif', 'w', driver='GTiff',
                                    height = array_to_export.shape[0], width = array_to_export.shape[1],
                                    count=1, dtype=str(perc.dtype),
                                    crs='EPSG:4326',
                                transform = transform)

        new_dataset.write(array_to_export, 1)
        new_dataset.close()

def agreement():
    """
    Subtracts historical percentiles from future to get change and sums agreement.
    Returns one raster with values from 0 to 10, representing the number of models
    that agree the percentile increases.
    """
    historical = rioxarray.open_rasterio('starting_data/percentile_historical.tif')
    
    CCSM4 = xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[0] +'.tif'), name = models[0])
    CESM1_BGC =  xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[1] +'.tif'), name = models[1])
    CMCC_CMS =  xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[2] +'.tif'), name = models[2])
    CanESM2 =  xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[3] +'.tif'), name = models[3])
    GFDL_CM3 =  xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[4] +'.tif'), name = models[4])
    HadGEM2_CC =  xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[5] +'.tif'), name = models[5])
    MIROC5 =  xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[6] +'.tif'), name = models[6])
    ACCESS1_0 =  xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[7] +'.tif'), name = models[7])
    HadGEM2_ES =  xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[8] +'.tif'), name = models[8])
    CNRM_CM5 =  xarray.DataArray.to_dataset(rioxarray.open_rasterio('tifs/'+scenario+'/percentile_' + models[9] +'.tif'), name = models[9])

    percentiles_rcp85 = xarray.merge(objects = [CCSM4, CESM1_BGC, CMCC_CMS, CanESM2, GFDL_CM3, HadGEM2_CC, MIROC5, ACCESS1_0, HadGEM2_ES, CNRM_CM5], compat = 'identical')

    # Subtract
    diff = percentiles_rcp85 - historical
    diff = diff.to_array()
    diff_trans = xarray.where(diff > 0, 1, 0)
    diff_sum = diff_trans.sum(dim = 'variable')

    diff_sum.plot()
    plt.show()

    # Export
    diff_sum.rio.to_raster('agreement_' + scenario + '.tif', dtype=np.int32)
    print('Exported tif for', scenario)


# Run all functions
get_percentile_array()
write_arrays_to_tifs()
agreement()