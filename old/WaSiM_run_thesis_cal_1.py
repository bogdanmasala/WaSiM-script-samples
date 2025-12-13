#   Script that automates control file modification and simulation in WaSiM,
#   with some output data processing steps
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

#set path to control folder
path_control = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Control'

#############################################################################################
def fun_substitute_line(oldline, newline, txtfile):   
    
    f1 = open(txtfile, 'r')
    f2 = open(txtfile + '.temp', 'w')
    
    try:
        for line in f1:
            
            if oldline in line:
                # You need to include a newline if you're replacing the whole line
                line = newline              
                f2.write(line + "\n")
            else:
                f2.write(line)

    except:
        print ("Line substitution in file: " + txtfile + " failed!!")
        
    f1.close()
    f2.close()
    
    try:    
        sh.copyfile(txtfile + '.temp',txtfile)
    except:
        print ("KOPPIERROR")
##################################################################################################
def fun_change_parameters(string, parameter):     

    param = parameter

    fun_substitute_line(string, string + '	' + str(param), path_control + '\\25_Martell_large_daily_20230523_1985_cal_1997_1.bash' )
    
    return param
###################################################################################################
def fun_change_parameters_rad(string, parameter):     

    param = parameter

    fun_substitute_line(string, '+' + str(param) + string, path_control + '\\25_Martell_large_daily_20230523_1985_cal_1997_1.bash' )
    
    return param
###################################################################################################

#load input matrix X from text file
X = np.loadtxt(os.path.join(path_control,'thesis_cal_X_round1.txt'))
#reduce decimal places to max. 2
np.round(X[:,0:6], decimals = 2, out = X[:,0:6])
np.round(X[:,6:10], decimals = 4, out = X[:,6:10])

#calibration parameters: copy text from each line for parameter of interest in the control file and define string variables
Ttrans = '$set $Ttrans		='
TRS = '$set $TRS			='
lwincorr = '$set $LWINcorr 		='
lwoutcorr = '$set $LWOUTcorr		='
MF = '$set $MF			='
VA_scal = '$set $VA_Scal		='
ice_min = '  	#0.0001			# radiation coefficient for ice_min  (for method 2)'
ice_max = '   					# radiation coefficient for ice_max  (for method 2)'
snow_min = ' 	#0.0001			# radiation coefficient for snow_min (for method 2)'
snow_max = '  	#0.00055 		# radiation coefficient for snow_max (for method 2)'

#define array with parameter names
X_parlabels_nonrad = [Ttrans, TRS, lwincorr, lwoutcorr, MF, VA_scal] 
X_parlabels_rad = [ice_min, ice_max, snow_min, snow_max]
#define number of parameters
num_parameters_nonrad = len(X_parlabels_nonrad)
#num_parameters_rad = len(X_parlabels_rad)

#define number of model runs/simulations and set equal to number of rows in input matrix X
num_modelruns = len(X[:,1])

#define names of WaSiM executables and control file
wasim = 'wasimvzo_100505_64.exe'
controlfile = '25_Martell_large_daily_20230523_1985_cal_1997_1.bash'

#get path to observations folder
path_obs = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Observation'

#create list of observed snow maps files and get dates
modis_list = []
modis_dir = os.path.join(path_obs, 'Snow_maps_Martelltal')
for modis_file in os.listdir(modis_dir):
    modis_list.append(os.path.join(modis_dir, modis_file))
    
#set Dataframe
snow_maps_obs = pd.DataFrame({'date' : [], 'filename' : []})    
snow_maps_obs['filename'] = modis_list

#date pattern search
date_pattern_modis = r'(\d{4}\-\d{2}\-\d{2})'
date_search_modis = []
for modis_file in modis_list:
    date_search_modis.append(re.search(date_pattern_modis,modis_file).group())
snow_maps_obs['date'] = pd.to_datetime(date_search_modis,format='mixed')
snow_maps_obs = snow_maps_obs.loc[snow_maps_obs['date']<'2015-10-01']

#for merging glacier extents with ssto raster, load 2006
path_init = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Init'
#2000 glacier cover
glc_2000_ds = rasterio.open(os.path.join(path_init,'glc_2000.asc'))
glc_2000 = glc_2000_ds.read(1)
glc_2000[glc_2000 == -9999] = 0
glc_2000[glc_2000 == 1] = 10
#2006 glacier cover
glc_2006_ds = rasterio.open(os.path.join(path_init,'glc_2006_larger_domain.asc'))
glc_2006 = glc_2006_ds.read(1)
glc_2006[glc_2006 == -9999] = 0
glc_2006[glc_2006 == 1] = 10
#2011 glacier cover
glc_2011_ds = rasterio.open(os.path.join(path_init,'glc_2011.asc'))
glc_2011 = glc_2011_ds.read(1)
glc_2011[glc_2011 == -9999] = 0
glc_2011[glc_2011 == 1] = 10

#open ezg raster for use as base in fgc calcs
path_input = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Input'
ezg_ds = rasterio.open(os.path.join(path_input,'dem_25_2006_large.ezg'))
ezg = ezg_ds.read(1)
ezg[ezg>0] = 0

#set mainpath, same as in control file, but without backslash
wasim_mainpath = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup'
    
#sims!
for runNumber in range(870,875): #possibility for running a part of the simulations
    #output folder and default output directory
    #define name of output folder
    output_folder_fname = 'output_thesis_cal_1997_round1_run_' + str(runNumber+1)
    #use path.join function to make new output folder path
    output_folder = os.path.join(wasim_mainpath, output_folder_fname)
    #change default output directory in control file
    oldline = '$set $DefaultOutputDirectory = $mainpath//'
    newline = oldline + output_folder_fname + '/'
    fun_substitute_line(oldline, newline, path_control + '\\' + controlfile)
    
    ###########################################################################
    #print error if output folders already exist, need delete folders function here!
    try:
        os.mkdir(output_folder)
    except OSError as error:
        print(error)
    #print runNumber+1 for easy reference in Spyder    
    print('run ' + str(runNumber+1))
    
    #change calibration parameters in control file
    for i in range(num_parameters_nonrad):
        #define variables for fun_change_parameters function
        string = X_parlabels_nonrad[i]
        parameter = (X[runNumber,i])
        fun_change_parameters((string), (parameter)) 
    for i in range(6,10):
        string = X_parlabels_rad[i-6]
        parameter = (X[runNumber,i])        
        fun_change_parameters_rad((string), (parameter))
            
    #copy modified control file of each parameter combination/simulation into corresponding output folder
    sh.copy(controlfile, wasim_mainpath + '\\' + output_folder_fname)
        
    #run WaSiM for first time period
    p = subprocess.Popen([wasim, controlfile])
    p.communicate()
    
    ###########################################################################
    #delete unneeded files
    for file_del in os.listdir(output_folder):
        if not file_del.startswith(('25_Martell_','glc','glmb','MA','ssto','special_','q')):
            os.unlink(os.path.join(output_folder,file_del))
        if file_del.startswith(('qu__stack','qi_','qifl')):
            os.unlink(os.path.join(output_folder,file_del))
    
    ###########################################################################
    #delete simulated ssto raster files that do not have matching dates with observations
    #create list of simulated ssto raster files
    sstodem_sim_list = []
    for output_file in os.listdir(output_folder):
        if output_file.startswith('sstodem_25_2006_large_'):
            sstodem_sim_list.append(os.path.join(output_folder, output_file))
    
    #create dataframe with simulated ssto raster dates and filenames
    ssto_sim = pd.DataFrame({'date' : [], 'filename' : []})    
    ssto_sim['filename'] = sstodem_sim_list  

    #collect all dates of simulated ssto rasters
    date_pattern_ssto = r'(\d{4}\.\d{2})'
    date_search_ssto = []
    for sstodem_file in sstodem_sim_list:
        date_search_ssto.append(re.search(date_pattern_ssto,sstodem_file).group())
    ssto_sim['date'] = pd.to_datetime(date_search_ssto,format='%y%m.%d')
    
    #delete ssto sim files for dates without observations
    ssto_delete = ~ssto_sim['date'].isin(snow_maps_obs['date'])
    #loop!
    for x in range(len(ssto_delete)):
        if ssto_delete.loc[x] == True:
            os.unlink(ssto_sim.loc[x, 'filename'])
            
    ###################### calculate fsc ######################################        
    #get new ssto raster list
    sstodem_new_list = []
    for output_file in os.listdir(output_folder):
        if output_file.startswith('sstodem_25_2006_large_'):
            sstodem_new_list.append(os.path.join(output_folder, output_file))
            
    #create dataframe of fsc values
    fsc_values = pd.DataFrame({'date': [], 'fsc_mean': [], 'fsc_sd': []})
    
    #create dataframe with left over simulated ssto raster dates and filenames
    ssto_sim_new = pd.DataFrame({'date' : [], 'filename' : []})    
    ssto_sim_new['filename'] = sstodem_new_list 
    #date search to fill date column
    date_pattern_ssto_new = r'(\d{4}\.\d{2})'
    date_search_ssto_new = []
    for sstodemnew_file in sstodem_new_list:
        date_search_ssto_new.append(re.search(date_pattern_ssto_new,sstodemnew_file).group())
    ssto_sim_new['date'] = pd.to_datetime(date_search_ssto_new,format='%y%m.%d')
    #plug dates in to fsc values table
    fsc_values['date'] = ssto_sim_new['date']
    
    #loop for each available simulated ssto raster
    for a in range(len(ssto_sim_new)):
        ssto_raster_ds = rasterio.open(ssto_sim_new.loc[a,'filename'])
        ssto_raster = ssto_raster_ds.read(1)
        #use different glcs for different years, also test with just using 2006 or so
        if ssto_sim_new.loc[a,'date'].year>2000 and ssto_sim_new.loc[a,'date'].year<2006:
            ssto_raster_glc = ssto_raster + glc_2000
        elif ssto_sim_new.loc[a,'date'].year>=2006 and ssto_sim_new.loc[a,'date'].year<2011:
            ssto_raster_glc = ssto_raster + glc_2006
        else:
            ssto_raster_glc = ssto_raster + glc_2011            
        #apply thresholds for snow detection
        ssto_raster_glc_final = ssto_raster_glc.copy()
        ssto_raster_glc_final[(ssto_raster_glc_final>=0) & (ssto_raster_glc_final<5)] = 0
        ssto_raster_glc_final[(ssto_raster_glc_final>=5) & (ssto_raster_glc_final<1000)] = 1
        ssto_raster_glc_final[(ssto_raster_glc_final>=1000) & (ssto_raster_glc_final<50000)] = 1
        
        ssto_raster_glc_final_calc = ssto_raster_glc_final.copy()
        ssto_raster_glc_final_calc[ssto_raster_glc_final_calc == -9999] = 'nan'

        fsc_values.iloc[a,1] = round(np.nanmean(ssto_raster_glc_final_calc),2)
        fsc_values.iloc[a,2] = round(np.nanstd(ssto_raster_glc_final_calc),2)

        #close ssto raster
        ssto_raster_ds.close()
        #print fsc values table
    fsc_values_table_path = output_folder + '\\table_fsc.csv'
    fsc_values.to_csv(fsc_values_table_path, sep = ',', index = False)
    
    #delete ssto rasters
    for file_del in os.listdir(output_folder):
        if file_del.startswith(('sstodem')) and file_del.endswith('m'):
            os.unlink(os.path.join(output_folder,file_del))

    ###########################################################################
    #glacier mass balances for Langenferner
    #delete simulated glmb files that do not have matching dates with observations
    glmb_sim_list = []
    for output_file in os.listdir(output_folder):
        if output_file.startswith('glmb') and output_file.endswith('2400'):
            glmb_sim_list.append(os.path.join(output_folder, output_file))
            
    #create dataframe with simulated glmb raster dates and filenames
    glmb_sim = pd.DataFrame({'date' : [], 'filename' : []})    
    glmb_sim['filename'] = glmb_sim_list             
            
    #collect all dates of simulated glmb rasters
    date_pattern_glmb = r'(\d{8})'
    date_search_glmb = []
    for glmb_file in glmb_sim_list:
        date_search_glmb.append(re.search(date_pattern_glmb,glmb_file).group())
    glmb_sim['date'] = pd.to_datetime(date_search_glmb,format='%Y%m%d')    
    #delete glmb files outside of Sep1/Oct1 or May balance measurement dates
    glmb_delete = ~glmb_sim['date'].dt.month.isin([4,5,9])
    #loop
    for y in range(len(glmb_delete)):
        if glmb_delete.loc[y] == True:
            os.unlink(glmb_sim.loc[y, 'filename'])    

    #check if rasters have to be kept or if only the fsc tables can be saved
    #ask about spatial comparison of obs and sim snow map rasters
    #save fsc tables for each snow threshold
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
        gc_raster = glc_raster + ezg
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

############################# anything else? ##################################
ezg_ds.close()
glc_2000_ds.close()
glc_2006_ds.close()        
glc_2011_ds.close()        
            
            
      
    
   



