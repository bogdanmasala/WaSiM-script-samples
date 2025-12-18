'''
This script adjusts a generated parameter set matrix where the snow and ice melt radiation coefficients are not
sampled with dependency and their ranges overlap. The parameter combinations where the coefficients for snow are
greater than the coefficients for ice are removed from the parameter set.
'''

#import libraries
import numpy as np
import os

#get and set path of working directory
path = os.getcwd() 
os.chdir(working directory)

#set directory paths and filepaths
parameter_set_matrix_filepath = 'parameter matrix filepath'
WaSiM_input_folder = 'input folder'

##############################################################################################################################################################################
def parameter_set_adjustment(parameter_set_matrix_filepath, input_folder): #substitute input folder for location of parameter matrix text files

    #load parameter set matrix from text file, with number of rows equal to number of parameter combinations
    parameter_set = np.loadtxt(parameter_set_matrix_filepath)
    #conditions: ice_min>snow_min and ice_max>snow_max
    parameter_set_adjusted = parameter_set.copy()

    #four coefficients for snow and ice radiation: snow_min, snow_max, ice_min, ice_max
    for a in range(len(parameter_set_adjusted)):
        if parameter_set[a, Ncolumn_snow_min] > parameter_set[a, Ncolumn_ice_min] or parameter_set[a, Ncolumn_snow_max] > parameter_set[a, Ncolumn_ice_max]:
            parameter_set_adjusted[a, :] = 0

    #portion of initially generated parameter combinations to be kept
    parameter_set_count = np.zeros(len(parameter_set))

    for b in range(len(parameter_set)):
        if parameter_set_adjusted[b, 1] != 0: #where parameter set values are greater than 0 or do not equal 0
            parameter_set_count[b] = 1
        else
            parameter_set_count[b] = 0

    #fraction of parameter combinations kept        
    fraction_kept = np.mean(parameter_set_count)

    #create array/matrix with only positive parameter combinations
    parameter_set_condition = np.all(parameter_set_adjusted != 0, 1) #where the adjusted parameter matrix does not have 0 values
    parameter_set_use = parameter_set[parameter_set_condition]

    #reduce decimal places to max. 2, reduce to N for snow and ice radiation coefficients
    np.round(parameter_set[:, 0:N], decimals = N_1, out = parameter_set[:, 0:N])
    np.round(parameter_set[:, Ncoefficient_columns:Ncoefficient_columns], decimals = N_2, out = parameter_set[:, Ncoefficient_columns:Ncoefficient_columns])

    #write adjusted parameter set matrix
    np.savetxt(os.path.join(input_folder, 'adjusted parameter set filename'), parameter_set_adjusted)      

    #return fraction of kept parameter combinations
    return fraction_kept
    
    #end function
##############################################################################################################################################################################

#run function
fraction_parameter_sets_kept = parameter_set_adjustment(parameter_set_matrix_filepath, input_folder)

#end

