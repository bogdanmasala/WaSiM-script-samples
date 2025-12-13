#   Script that automates some WaSiM simulation tasks for a 
#   glacier initialization in 1850
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
path_input = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup\\Input'
wasim_mainpath = 'C:\\Users\\Asus\\Documents\\Thesis\\WaSiM_setup'

######################## functions for changing lines in control files ########
###############################################################################
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
###############################################################################
def fun_change_parameters(string, parameter):     

    param = parameter

    fun_substitute_line(string, string + '	' + str(param), path_control + '\\25_Martell_large_daily_20230523_1850_2.bash' )
    
    return param
###############################################################################
def fun_change_parameters_rad(string, parameter):     

    param = parameter

    fun_substitute_line(string, '+' + str(param) + string, path_control + '\\25_Martell_large_daily_20230523_1850_2.bash' )
    
    return param
###############################################################################
######################### control file parameters #############################
#load input matrix X from text file
X = np.loadtxt(os.path.join(path_control,'thesis_X_1850_use.txt'))
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

###############################################################################
################## set up zeroed ezg raster for fgc and glmb  #################
#open ezg raster for use as base in fgc calcs
fname_ezg = 'dem_25_2006_large.ezg'
ezg_ds = rasterio.open(os.path.join(path_input,fname_ezg))
ezg = ezg_ds.read(1)
ezg_0 = ezg.copy()
ezg_0[ezg_0 >= 0] = 0

###############################################################################
#define names of WaSiM executables and control file
wasim = 'wasimvzo_100505_64.exe'
controlfile = '25_Martell_large_daily_20230523_1850_2.bash'

################################### sims! #####################################
for runNumber in range(15,16): #run sims in batches
    ############ create output folders and modify control file ################
    #output folder and default output directory
    #define name of output folder
    output_folder_fname = 'output_thesis_1850_run_' + str(runNumber+1)
    #use path.join function to make new output folder path
    output_folder = os.path.join(wasim_mainpath, output_folder_fname)
    #change default output directory in control file
    oldline = '$set $DefaultOutputDirectory = $mainpath//'
    newline = oldline + output_folder_fname + '/'
    fun_substitute_line(oldline, newline, path_control + '\\' + controlfile)
    
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
    
    ###################### delete unneeded files ##############################
    for file_del in os.listdir(output_folder):
        if not file_del.startswith(('25_Martell_','glc','MA','ssto','special_','q')):
            os.unlink(os.path.join(output_folder,file_del))
        if file_del.startswith(('qu__stack','qi_','qifl')):
            os.unlink(os.path.join(output_folder,file_del))

    ###########################################################################
    ########################### fgc calculations ##############################
    #get list of glc rasters
    glc_sim_list = []
    for output_file in os.listdir(output_folder):
        if output_file.endswith('2400') and output_file.startswith('glc_'):
            glc_sim_list.append(os.path.join(output_folder, output_file))
    
    #initialize dataframe for date and fgc mean 
    fgc_values = pd.DataFrame({'date': [],'fgc': []})
    
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
        fgc_values.loc[b,'fgc'] = round(np.nanmean(gc_raster),2)
        glc_raster_ds.close()
    
    #write table of fgc values
    fgc_values_table_path = output_folder + '\\table_fgc.csv'
    fgc_values.to_csv(fgc_values_table_path, sep = ',', index = False)     
    
    #delete glc and glid files between 1850 and 1997, not inclusive
    #create dataframe with simulated glc raster dates and filenames
    glc_sim = pd.DataFrame({'date' : [], 'filename' : []})    
    glc_sim['filename'] = glc_sim_list             
    #collect all dates of simulated glc rasters
    glc_sim['date'] = fgc_values['date']    
    #delete with loop
    for x in range(len(glc_sim)):
        if glc_sim.loc[x,'date'].year>1850 and glc_sim.loc[x,'date'].year<1997:
            os.unlink(glc_sim.loc[x, 'filename'])  
    ###########################################################################
    
#close all open rasters
ezg_ds.close()

############################# anything else? ##################################



