#   Script for comparing the number and timing of snow days from different model outputs
#   with that of an observational series
import spotpy
import numpy as np
import os
import pandas as pd
import subprocess
from datetime import datetime as dt
import re
import math
path = os.getcwd() 
import sys
sys.path.insert(0, path)
import shutil as sh
from dateutil.parser import parse
import matplotlib.pyplot as plt   

#functions and for loops need to be added for efficiency

############################################################
#calculate days of increasing SWE in observations
######################################################################
obs_file = '11200BM_StauseeZufrittBeobachter_DigadiGioverettoOsservatore_HS_SWE.Beobachter_0900Day.Cmd.csv'
obs_dir = 'C:\\Users\\Asus\\Documents\\StudyProject2023S\\WaSiM_setuptest\\Observation\\Snow'
obs_path = os.path.join(obs_dir, obs_file)
zufritt_obs = pd.read_csv(obs_path, sep = ',', engine = 'python', parse_dates=True)
O = zufritt_obs[9488:10949]

O_snowdays = np.zeros(len(O))
a=0
for a in range(len(O_snowdays)-1):
    if float(O.loc[a+9489,'SWE_Pistocchi'])>float(O.loc[a+9488,'SWE_Pistocchi']):
        O_snowdays[a+1] = 1
    else:
        O_snowdays[a+1] = 0
  
############################################################
#use sday output, using subcatchment 1
######################################################################
sday_dir = 'C:\\Users\\Asus\\Documents\\StudyProject2023S\\WaSiM_setuptest\\Test_output'
sday_file = 'sdaydem_25_2006_large.stat'
sday_path = os.path.join(sday_dir,sday_file)
sday = np.loadtxt(sday_path,dtype=str,skiprows=2,usecols=(0,1,2,3,4))

sday_snowdays = np.zeros(len(O))
c=0
for c in range(len(O)-1):
    if sday[c+2,4]>sday[c+1,4]:
        sday_snowdays[c+1] = 1
    else:
        sday_snowdays[c+1] = 0

############################################################
#use snow rate output, using subcatchment 1
######################################################################
snow_dir = 'C:\\Users\\Asus\\Documents\\StudyProject2023S\\WaSiM_setuptest\\Test_output'
snow_file = 'snowdem_25_2006_large.stat'
snow_path = os.path.join(snow_dir,snow_file)
snow = np.loadtxt(snow_path,dtype=str,skiprows=2,usecols=(0,1,2,3,4))

snow_snowdays = np.zeros(len(O))
d=0
for d in range(len(O)-1):
    if float(snow[d+1,4])>0:
        snow_snowdays[d] = 1
    else:
        snow_snowdays[d] = 0

############################################################
#create snow rate special output at zufritt station?
######################################################################

############################################################
#fsc timing?
######################################################################
 
############################################################
#calculate days of increasing SWE in 'otherSnowData' output, using  zufritt station location
###################################################################### 
ssto_dir = 'C:\\Users\\Asus\\Documents\\StudyProject2023S\\WaSiM_setuptest\\output_run_glacierSA_55'
ssto_file = 'special_output_otherSnowData.stat2011'  
ssto_path = os.path.join(ssto_dir,ssto_file)
ssto = np.loadtxt(ssto_path,dtype=str,skiprows=2,usecols=(0,1,2,3,4,6,8))

ssto_snowdays = np.zeros(len(O))
b=0
for b in range(len(O)-1):
    if float(ssto[b+2,6])>float(ssto[b+1,6]):
        ssto_snowdays[b+1] = 1
    else:
        ssto_snowdays[b+1] = 0

############################################################
#station 2401.09 from WRF data
######################################################################    
WRF2401_dir = 'C:\\Users\\Asus\\Documents\\StudyProject2023S\\WaSiM_setuptest\\Input'  
WRF2401_file = 'WRF_Martelltal_daily_precipitation_1850_2015.txt'
WRF2401_path = os.path.join(WRF2401_dir,WRF2401_file)
WRF2401 = np.loadtxt(WRF2401_path,dtype=str,skiprows=1,usecols=(0,1,2,3,23))
WRF2401 = WRF2401[57377:58838]

WRF2401_precdays = np.zeros(len(O))
f=0
for f in range(len(O)-1):
    if float(WRF2401[f+1,4])>float(WRF2401[f,4]):
        WRF2401_precdays[f+1] = 1
    else:
        WRF2401_precdays[f+1] = 0    

#set approximate summer non-snow precipitation to 0 using ssto numbers, don't take it too seriously
#true timeline
WRF2401_precdays_adj = WRF2401_precdays
WRF2401_precdays_adj[247:332] = 0
WRF2401_precdays_adj[618:731] = 0
WRF2401_precdays_adj[1004:1096] = 0
WRF2401_precdays_adj[1370:1432] = 0

############################################################
#make summary dataframe
######################################################################
snowdays = pd.DataFrame({'date':[],'obs':[],'ssto':[],'sday':[],'snow':[]})
snowdays['date'] = pd.date_range(start='2006-10-01',end='2010-09-30',freq ='D')
snowdays['date'] = snowdays['date'].dt.date
snowdays['obs'] = O_snowdays
snowdays['ssto'] = ssto_snowdays
snowdays['sday'] = sday_snowdays
snowdays['snow'] = snow_snowdays

############################################################
#RMSE for timing between obs and ssto
######################################################################
#residual
residual_ssto = np.zeros((len(snowdays),1))
e=0
for e in range(len(residual_ssto)):
    residual_ssto[e,0] = snowdays.loc[e,'obs'] - snowdays.loc[e,'ssto']
residual_ssto_sq = np.square(residual_ssto)    
   
#RMSE (root mean square error)
RMSE_ssto = math.sqrt((residual_ssto_sq).mean())

############################################################
#RMSE for timing between obs and sday
######################################################################
#residual
residual_sday = np.zeros((len(snowdays),1))
e=0
for e in range(len(residual_sday)):
    residual_sday[e,0] = snowdays.loc[e,'obs'] - snowdays.loc[e,'sday']
residual_sday_sq = np.square(residual_sday)    
   
#RMSE (root mean square error)
RMSE_sday = math.sqrt((residual_sday_sq).mean())

############################################################
#RMSE for timing between obs and snow
######################################################################
#residual
residual_snow = np.zeros((len(snowdays),1))
e=0
for e in range(len(residual_snow)):
    residual_snow[e,0] = snowdays.loc[e,'obs'] - snowdays.loc[e,'snow']
residual_snow_sq = np.square(residual_snow)    
   
#RMSE (root mean square error)
RMSE_snow = math.sqrt((residual_snow_sq).mean())

############################################################
#RMSE for timing between obs and WRF2401
######################################################################
#residual
residual_WRF2401 = np.zeros((len(snowdays),1))
e=0
for e in range(len(residual_WRF2401)):
    residual_WRF2401[e,0] = snowdays.loc[e,'obs'] - WRF2401_precdays_adj[e]
residual_WRF2401_sq = np.square(residual_WRF2401)    
   
#RMSE (root mean square error)
RMSE_WRF2401 = math.sqrt((residual_WRF2401_sq).mean())

############################################################
#make some plots for presentations
######################################################################
plt.plot(snowdays['date'], snowdays['ssto'], 'r.', markersize=1)
plt.plot(snowdays['date'], snowdays['obs'], 'b.', markersize=1)
plt.plot(snowdays['date'], snowdays['snow'], 'g.', markersize=1)
plt.plot(snowdays['date'], snowdays['sday'], 'y.', markersize=1)








