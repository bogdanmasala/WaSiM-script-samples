'''
This script conducts goodness-of-fit testing for snow measurement station data when model-generated weather data
was used as input to the WaSiM model.
'''

#import libraries
import numpy as np
import pandas as pd
import os
import re
from datetime import datetime as dt
from scipy.stats import kstest
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt   

#get and set path of working directory
path = os.getcwd() 
os.chdir(working directory)

#set directory paths and filenames
WaSiM_observations_folder = 'observations_folder'
output_folder = 'output'

##############################################################################################################################################################################
def snow_station_ecdf(observations_folder, snow_station_observations_file, snow_station_output_file, output_folder):
    
    #open observations time series
    observations = pd.read_csv(os.path.join(observations_folder, snow_station_observations_file), sep = ';', engine = 'python', parse_dates = True)
    
    #take only the observations that correspond to the simulation period
    observations_simulation_period = observations[Nstart, Nend]
    
    #initialize dataframe for comparing observations and simulations
    comparison = pd.DataFrame({'date': [], 'observations': [], 'output': []})
    #date_pattern = r'(\d{8}\)'
    comparison['date'] = pd.to_datetime(observations_simulation_period['NameofDateColumn'], format = '%Y%m%d')
    comparison['date'] = pd.to_datetime(comparison['date']).dt.date
    comparison['observations'] = observations_simulation_period['NameofObservationsColumn']
    
    #open output data, or simulations, file
    simulations = np.loadtxt(os.path.join(output_folder, snow_station_output_file), dtype = str, skiprows = Nskiprows, usecols = (N1...Nn))
    
    #additional processing needed if observational data does not include all days within period
    #create dataframe to assign dates to simulation observations
    simulations_dates = pd.DataFrame({'date': [], 'output': []})
    simulations_dates['date'] = pd.date_range(start = 'YYY-MM-DD', end = 'YYY-MM-DD', freq ='D')
    simulations_dates['date'] = simulations_dates['date'].dt.date
    simulations_dates['output'] = simulations[Nstart:Nend, Ncolumn]
    #more efficient process needed below
    for x in range(len(comparison)):
        for y in range(len(simulations_dates)):
            if comparison.loc[x,'date'] == simulations_dates.loc[y,'date']:
                comparison.loc[x,'output'] = float(simulations_dates.loc[y,'output'])
    
    #K-S test
    KStest = np.array([kstest(comparison['output'], comparison['observations']), ks_2samp(comparison['output'], comparison['observations'])]) 
   
    #optional step: plot time series
    SWE_compare = SWE_compare.set_index('date')

    #optional step: plot ecdf's
    fig, ax = plt.subplots(figsize = (6, 4))
    ax.set_title('ECDF')
    ax.set_ylabel('fraction non-exceedance')
    ax.set_xlabel('SWE (mm)')
    ax.axis([0, max(comparison['output']), 0.0, 1.0])
    ax.grid(True)
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    for column in ['output', 'observations']
        ax.ecdf(comparison[column], label = 'SWE')    
    ax.legend()
   
    #return K-S test results
    return KStest
    
    #end function
##############################################################################################################################################################################    
    
#run function
KStest = def snow_station_ecdf(observations_folder, snow_station_observations_file, snow_station_output_file, output_folder)


#end




