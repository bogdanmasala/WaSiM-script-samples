#   Script for computing fractional glacier-covered area of your catchment
import spotpy
import numpy as np
import os
import pandas as pd
import subprocess
from datetime import datetime
path = os.getcwd() 
import sys
sys.path.insert(0, path)
import shutil as sh
import re
import rasterio

######################### set paths to needed folders #########################
path_control = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Control'
path_obs = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Observation'
path_init = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Init'
path_langerner = os.path.join(path_init,'langenferner') 
path_input = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Input'
wasim_mainpath = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\supplement runs'
#change working directory to control folder
os.chdir(path_control)

###############################################################################
################## set up zeroed ezg raster for fgc #################
#open ezg raster for use as base in fgc calcs
fname_ezg = 'dem_25_2006_large.ezg'
ezg_ds = rasterio.open(os.path.join(path_input,fname_ezg))
ezg = ezg_ds.read(1)
ezg_0 = ezg.copy()
ezg_0[ezg_0 >= 0] = 0


output_folder_fname = 'output_run_thesis_WRF_firn_mod'
#use path.join function to make new output folder path
output_folder = os.path.join(wasim_mainpath, output_folder_fname)
#change default output directory in control file


###########################################################################
########################### fgc calculations ##############################
#get list of glc rasters
glc_sim_list = []
for output_file in os.listdir(output_folder):
    if output_file.endswith('2400') and output_file.startswith('glc_'):
        glc_sim_list.append(os.path.join(output_folder, output_file))

#initialize dataframe for date and fgc mean 
fgc_values = pd.DataFrame({'date': [],'fgc_mean': []})

#collect dates of simulated glc rasters
date_pattern_glc = r'(\d{8})'
date_search_glc = []
for glc_file in glc_sim_list:
    date_search_glc.append(re.search(date_pattern_glc,glc_file).group())
fgc_values['date'] = pd.to_datetime(date_search_glc,format='%Y%m%d')

#calculate fgc means and fill in fgc_values table
for b in range(len(glc_sim_list)):
    glc_raster_ds = rasterio.open(glc_sim_list[b])
    glc_raster = glc_raster_ds.read(1)
    glc_raster[glc_raster<0] = 0
    
    #add glc raster to base raster
    gc_raster = glc_raster + ezg_0
    gc_raster[gc_raster<0] = 'nan'
    
    #take NaN mean of glacier cover raster
    fgc_values.loc[b,'fgc_mean'] = round(np.nanmean(gc_raster),2)
    glc_raster_ds.close()

#write table of fgc values
fgc_values_table_path = output_folder + '\\table_fgc.csv'
fgc_values.to_csv(fgc_values_table_path, sep = ',', index = False)     

#don't delete glc files! or maybe do? dont need spatial fit for first sample runs
#delete fgc rasters after all
for file_del in os.listdir(output_folder):
    if file_del.startswith(('glc_')) and file_del.endswith('2400'):
        os.unlink(os.path.join(output_folder,file_del))  

ezg_ds.close()

