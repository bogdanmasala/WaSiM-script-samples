'''
This script creates binary snow cover rasters, using the snow depth rasters of the WaSiM model output. It then computes 
fractional snow-covered area for the catchment, as well as standard deviation of binary snow cover values.
'''

#import libraries
import numpy as np
import pandas as pd
import rasterio
import os
import re
from datetime import datetime as dt

#get and set path of working directory
path = os.getcwd() 
os.chdir(working directory)

#set directory paths and filenames
WaSiM_input_folder = 'input_folder'
glacier_base_raster = 'glacier_base.asc'
output_folder = 'output_folder'

##############################################################################################################################################################################
def fractional_snow_covered_area(input_folder, glacier_base_raster, output_folder):
    
    #Read a glacier extents raster from an appropriate year during the simulation period, glacierized cells should be counted as snow-covered cells
    #assign an SWE threshold value to count glacierized cells as snow-covered, i.e. change 1's to 10's
    glacier_base_ds = rasterio.open(os.path.join(input_folder, glacier_base_raster))
    glacier_base = glacier_base_ds.read(1)
    #set no-data values to 0
    glacier_base[glacier_base == -9999] = 0
    glacier_base[glacier_base == 1] = 10

    #make list of filepaths of ssto output files
    ssto_list = []
    for search_file in os.listdir(output_folder):
        if search_file.startswith('sstodem...'):
            ssto_list.append(os.path.join(output_folder, search_file))

    #initialize dataframe with date of output, fraction of catchment area covered in snow, and standard deviation of binary snow cover grid values
    fraction_snow_cover = pd.DataFrame({'date': [], 'fraction_snow_cover': [], 'standard_deviation_snow_cover': []})
    
    #collect dates of ssto files
    date_pattern = r'(\d{4}\.\d{2})'
    date_search = []
    for date_file in ssto_list:
        date_search.append(re.search(date_pattern, date_file).group())
    fraction_snow_cover['date'] = pd.to_datetime(date_search, format = '%y%m.%d')
    fraction_snow_cover['date'] = pd.to_datetime(fraction_snow_cover['date']).dt.date

    #compute fractional snow-covered area and standard deviation, and fill in table
    for a in range(len(ssto_list)):
    
        ssto_raster_ds = rasterio.open(ssto_list[a])
        ssto_raster = ssto_raster_ds.read(1)
        ssto_raster_with_glaciers = ssto_raster + glacier_base
        
        #compute binary snow cover rasters (1 = snow-covered cell, 0 = not snow-covered cell) using an SWE threshold of 5mm, or other
        ssto_raster_binary = ssto_raster_with_glaciers.copy()
        ssto_raster_binary[(ssto_raster_binary >= 0) & (ssto_raster_binary < XSWEthreshold)] = 0
        ssto_raster_binary[(ssto_raster_binary >= XSWEthreshold) & (ssto_raster_binary < 1000)] = 1
        ssto_raster_binary[(ssto_raster_binary >= 1000) & (ssto_raster_binary < 50000)] = 1

        #make copy for calculations
        ssto_raster_binary_calculation = ssto_raster_binary.copy()
        ssto_raster_binary_calculation[ssto_raster_binary_calculation == -9999] = 'nan'
        
        #compute NaN (not a number) mean and standard deviation of computed binary snow cover rasters
        fraction_snow_cover.iloc[x, 1] = round(np.nanmean(ssto_raster_binary_calculation), 2)
        fraction_snow_cover.iloc[x, 2] = round(np.nanstd(ssto_raster_binary_calculation), 2)

        #save profile of ssto output rasters
        kwds_ssto = ssto_raster_ds.profile

        #write binary snow map raster
        #define file path for generated raster
        snow_map_raster_path =  output_folder + '\\Snow_maps_' + str(fraction_snow_cover.loc[x,'date']) + '.asc'
        snow_map_raster_file = rasterio.open(snow_map_raster_path, mode = 'w', **kwds_ssto)
        snow_map_raster_file.write(ssto_raster_binary, 1)
        #close open raster file
        snow_map_raster_file.close()
    
    #write table of snow cover fraction and standard deviation
    fraction_snow_cover_table_path = output_folder + '\\table_snow_fractions.csv'
    fraction_snow_cover.to_csv(fraction_snow_cover_table_path, sep = ',', index = False)
        
    #close open raster file
    glacier_base_ds.close()

    #end function
##############################################################################################################################################################################
                
#run function
fractional_snow_covered_area(input_folder, glacier_base_raster, output_folder)

