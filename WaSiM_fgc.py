'''
This script reads the binary glacier cover files from an output folder of a WaSiM model run, and computes the fraction of 
the catchment area covered in glaciers, or the fractional glacier-covered area ("fgc"). A table summarizing the computations
is written into the output folder, and the binary glacier cover files can then be deleted to save space.
'''

#import libraries
import numpy as np
import pandas as pd
import rasterio
import os
import re

#get and set path of working directory
path = os.getcwd() 
os.chdir(working directory)

#set directory paths and filenames
WaSiM_control_folder = 'control_folder'
WaSiM_observations_folder = 'observations_folder'
WaSiM_initialization_folder = 'init_folder'
WaSiM_input_folder = 'input_folder'
WaSiM_main_directory = 'main_directory'
output_folder = 'output'

##############################################################################################################################################################################
def fractional_glacier_covered_area(input_folder, subcatchment_identification_file, output_folder):
    
    #set up zeroed subcatchments raster
    subcatchments_ds = rasterio.open(os.path.join(input_folder, subcatchment_identification_file))
    subcatchments = subcatchments_ds.read(1)
    subcatchments_zeroed = subcatchments.copy()
    subcatchments_zeroed[subcatchments_zeroed > 0] = 0
    
    #make list of filepaths of output glacier cover rasters
    glacier_cover_list = []
    for search_file in os.listdir(output_folder):
        if search_file.endswith('file ending') and search_file.startswith('file start'):
            glacier_cover_list.append(os.path.join(output_folder, search_file))
            
    #initialize dataframe with date of output and fraction of catchment area covered in glacier
    fraction_glacier_cover = pd.DataFrame({'date': [], 'fraction_glacier_cover': []})
    
    #collect dates of output glacier cover rasters
    date_pattern = r'(d{8})'
    date_search = []
    for date_file in glacier_cover_list:
        date_search.append(re.search(date_pattern, date_file).group())
    fraction_glacier_cover['date'] = pd.to_datetime(date_search, format='%Y%m%d')

    #compute fractional glacier-covered area and fill in table
    for a in range(len(glacier_cover_list)):
        glacier_cover_ds = rasterio.open(glacier_cover_list[a])
        glacier_cover = glacier_cover_ds.read(1)
        #set non-data values to 0
        glacier_cover[glacier_cover < 0] = 0
        
        #add zeroed subcatchments raster
        glacier_cover = glacier_cover + subcatchments_zeroed
        glacier_cover[glacier_cover < 0] = 'nan'
        
        #compute NaN (not a number) mean of glacier cover raster
        fraction_glacier_cover.loc[a, 'fraction_glacier_cover'] = round(np.nanmean(glacier_cover), 2)
        
        #close opened raster file
        glacier_cover_ds.close()
        
    #write table of fractional glacier-covered area values
    fraction_glacier_cover_table_path = output_folder + '\\table_glacier_fractions.csv'
    fraction_glacier_cover.to_csv(fraction_glacier_cover_table_path, sep = ',', index = False)     

    #can delete output glacier cover rasters
    #delete fgc rasters after all
    for delete_file in os.listdir(output_folder):
        if delete_file.startswith(('file ending')) and file_del.endswith('file start'):
            os.unlink(os.path.join(output_folder, delete_file))  

    #close open raster files
    subcatchments_ds.close()    

    #end function
##############################################################################################################################################################################

#run function
fractional_glacier_covered_area(input_folder, subcatchment_identification_file, output_folder)



