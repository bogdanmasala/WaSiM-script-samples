'''
Early and basic version of a script that runs a batch a WaSiM simulations. Needs an array of parameter combinations,
in my case, created using the SAFE Toolbox. Does not include data processing. The script has three basic steps:
1.  Creates an output folder for each simulation
2.  Modifies the WaSiM control file, creates a copy in each simulation folder, 
3.  Runs the WaSiM executable with the corresponding control file.

'''

#   First draft of a script that automates control file modification, 
#   output folder creation, and simulation runs in WaSiM
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

#############################################################################################
def fun_substitute_line(oldline, newline, txtfile):   
    
    f1 = open(txtfile, 'r')
    f2 = open(txtfile + '.temp', 'w')
    
    try:
        for line in f1:
            
            if oldline in line:
#                print oldline
                # You need to include a newline if you're replacing the whole line
                line = newline              
#                print newline
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
#def fun_change_parameters(gridName, runNumber,string, parameter):   
def fun_change_parameters(string, parameter):     

    #LOAD ALIAS
    #FOLDERS, FILES = fun_load_alias(gridName, runNumber) 
    
    param = parameter
    
    #CHANGE CONTROLFILE
    fun_substitute_line(string, string + '	' + str(param), path + '\\25_Martell_large_daily_20230523_2006_scripttest.bash' )
    
    return param
#########################################################################################

#load input matrix X from text file
X = np.loadtxt('')

#make test input matrix X
X = np.array([[0.1, 0.6, 0.8, 0.8, 0.8, 1, 2600, 26],[5, 1, 1.3, 1.3, 1.5, 3, 3000, 30]])

#copy text from each line for parameter of interest in the control file and define string variables
Ttrans = '$set $Ttrans		='
TRS = '$set $TRS			='
lwincorr = '$set $LWINcorr 		='
lwoutcorr = '$set $LWOUTcorr		='
spfactor = '$set $spfactor		='
MF = '$set $MF			='
ELA = '$set $ELA			='
VA_scal = '$set $VA_Scal		='
itera = '$set $iter			=' 

#define array with parameter names
X_parlabels = [Ttrans, TRS, lwincorr, lwoutcorr, spfactor, MF, ELA, VA_scal]
#define number of parameters
num_parameters = len(X_parlabels)

#define number of model runs/simulations and set equal to number of rows in input matrix X
num_modelruns = len(X[:,1])

#define names of WaSiM executables and control file
wasim = 'wasimvzo_100505_64.exe'
controlfile = '25_Martell_large_daily_20230523_2006_scripttest.bash'

#initialize run number and other loop
runNumber = 0
i = 0

for runNumber in range(num_modelruns):
    #set mainpath, same as in control file, but without backslash
    wasim_mainpath = 'C:\\Users\\Asus\\Documents\\StudyProject2023S\\WaSiM_setuptest'
    #define name of output folder for each simulation
    output_folder = 'output_run_' + str(runNumber+1)
    #use path.join function
    run_folder = os.path.join(wasim_mainpath, output_folder)
    
    #change default output directory in control file for each simulation
    oldline = '$set $DefaultOutputDirectory = $mainpath//'
    newline = oldline + output_folder + '/'
    fun_substitute_line(oldline, newline, path + '\\25_Martell_large_daily_20230523_2006_scripttest.bash')
    
    #print error if output folders already exist, need delete folders function here!
    try:
        os.mkdir(run_folder)
    except OSError as error:
        print(error)
        
    #optional: change simulation time in control file    
    
    for i in range(num_parameters):
        #define variables for fun_change_parameters function
        string = X_parlabels[i]
        parameter = (X[runNumber,i])
        fun_change_parameters((string), (parameter))
        
        #copy control file of each simulation into corresponding output folder
        sh.copy(controlfile, wasim_mainpath + '\\' + output_folder)
    
    #run WaSiM
    p = subprocess.Popen([wasim, controlfile])
    p.communicate()

