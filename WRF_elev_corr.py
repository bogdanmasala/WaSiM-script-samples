'''
Precipitation series from WRF sometimes need to be adjusted based on elevation. This script creates a 
precipitation correction raster using a formula estimated by a participant on the SEHAG project.

'''

import numpy as np
import os
import re
import rasterio
path = os.getcwd() 
import sys

#optimization with dictionaries and functions needed

#import elevation raster
path_input = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Input'
input_dhm_ds = rasterio.open(os.path.join(path_input,'dem_25_2006_large.dhm'))
input_dhm = input_dhm_ds.read(1)

########################### using dhm raster ##################################
#### not needed anymore, use only catchment area for elevation correction #####
#mean and median elevation, percentiles, standard deviation
#change -9999 to nan for  statistics calculations
input_dhm_calc = input_dhm.copy()
input_dhm_calc[input_dhm_calc == -9999] = 'nan'

elev_mean = np.nanmean(input_dhm_calc)
elev_median = np.nanmedian(input_dhm_calc)
elev_1pc = round(np.nanpercentile(input_dhm_calc, 1),2)
elev_10pc = round(np.nanpercentile(input_dhm_calc, 10),2)
elev_50pc = round(np.nanpercentile(input_dhm_calc, 50),2)
elev_90pc = round(np.nanpercentile(input_dhm_calc, 90),2)
elev_100pc = round(np.nanpercentile(input_dhm_calc, 100),2)
elev_std = round(np.nanstd(input_dhm_calc),2)

#regression equation
#50% precipitation change for 1000m change in elevation
#make copy of raster
input_dhm_corr = input_dhm.copy()
input_dhm_corr[input_dhm_corr == -9999] = 'nan'
input_dhm_corr[input_dhm_corr > 0] = 0
#calculate correction factors grid
input_dhm_corr[:,:] = 1 + 0.5*(input_dhm_calc[:,:] - elev_50pc)/(1000)

#statistics of correction factors grid
dhm_corr_mean = np.nanmean(input_dhm_corr)
dhm_corr_1pc = round(np.nanpercentile(input_dhm_corr, 1),2)
dhm_corr_10pc = round(np.nanpercentile(input_dhm_corr, 10),2)
dhm_corr_90pc = round(np.nanpercentile(input_dhm_corr, 90),2)
dhm_corr_100pc = round(np.nanpercentile(input_dhm_corr, 100),2)
dhm_corr_std = round(np.nanstd(input_dhm_corr),2)

############################## load ezg input grid ############################

#raster with catchment area
#load ezg grid
input_ezg_ds = rasterio.open(os.path.join(path_input,'dem_25_2006_large.ezg'))
input_ezg = input_ezg_ds.read(1)
#make copy for further calculations and set non-zero cells to 0
input_ezg_calc = input_ezg.copy()
#input_ezg_calc =input_ezg_calc.astype(dtype = float)
input_ezg_calc[input_ezg_calc>0]=0

########################### dhm limited to catchment area #####################

#create new elevation raster limited to catchment area, using raster ezg
dhm_ezg = input_dhm.copy()
dhm_ezg = dhm_ezg + input_ezg_calc
dhm_ezg[dhm_ezg < 0] = -9999

#calculate new median elevation
#make copy of input_dhm_ezg for statistics calculations
dhm_ezg_calc = dhm_ezg.copy()
dhm_ezg_calc[dhm_ezg_calc == -9999] = 'nan'
elev_ezg_50pc = round(np.nanpercentile(dhm_ezg_calc, 50),2)

#make copy of dhm_ezg raster for calculating correction factors
dhm_ezg_corr = dhm_ezg.copy()

#regression equation again
#50% precipitation change for 1000m change in elevation
#calculate correction factors grid
dhm_ezg_corr[:,:] = 1 + 0.5*(dhm_ezg[:,:] - elev_ezg_50pc)/(1000)
#set no data cells to -9999
dhm_ezg_corr[dhm_ezg_corr < 0] = -9999
#np.round doesn't change decimal places in written file
#dhm_ezg_corr = np.round(dhm_ezg_corr,decimals=4)

#statistics of correction factors grid
#make copy of dhm_ezg_corr for statistics calculations
dhm_ezg_corr_calc = dhm_ezg_corr.copy()
dhm_ezg_corr_calc[dhm_ezg_corr_calc == -9999] = 'nan'
dhm_ezg_corr_mean = np.nanmean(dhm_ezg_corr_calc)
dhm_ezg_corr_1pc = round(np.nanpercentile(dhm_ezg_corr_calc, 1),2)
dhm_ezg_corr_10pc = round(np.nanpercentile(dhm_ezg_corr_calc, 10),2)
dhm_ezg_corr_90pc = round(np.nanpercentile(dhm_ezg_corr_calc, 90),2)
dhm_ezg_corr_100pc = round(np.nanpercentile(dhm_ezg_corr_calc, 100),2)
dhm_ezg_corr_std = round(np.nanstd(dhm_ezg_corr_calc),2)

########################## for loop array calculation as check ################
# nah!
########################### profiles and writing rasters ######################

#save profile of dhm and ezg rasters
kwds_dhm = input_dhm_ds.profile
kwds_ezg = input_ezg_ds.profile

#write rasters
#needs to be float datatype, decimals missing on first try
#need grid for each month
month_extensions = np.array([1,2,3,4,5,6,7,8,9,10,11,12], dtype = str)
for a in range(9):
    #print(a)
    dhm_ezg_corr_file = rasterio.open(os.path.join(path_input,'WRF_elevation_correction.0' + month_extensions[a]), mode = 'w', **kwds_dhm)
    dhm_ezg_corr_file.write(dhm_ezg_corr,1)
    dhm_ezg_corr_file.close()
for a in range(9,12):
    #print(a)
    dhm_ezg_corr_file = rasterio.open(os.path.join(path_input,'WRF_elevation_correction.' + month_extensions[a]), mode = 'w', **kwds_dhm)
    dhm_ezg_corr_file.write(dhm_ezg_corr,1)
    dhm_ezg_corr_file.close()    

#close rasters
input_dhm_ds.close()
input_ezg_ds.close()

########## plot differences in qgko, SWE at zufritt and zufallhuette ##########

















