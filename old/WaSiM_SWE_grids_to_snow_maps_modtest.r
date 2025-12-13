#	R script that converts of simulated WaSiM SWE grids to binary snow maps

rm(list = ls())
# load required packages
library(lubridate)
library(zoo)
library(raster)
library(stringr)
library(rgdal)
library(matrixStats)
library(data.table)

# Case with one calibration variable
setwd("C:/Users/Asus/Documents/StudyProject2023S/testoutputfolders")

# Snow-rain threshold parameters
output_folders = seq(1,6,1) #sequence with number of output folders

# Folder names and parameter configuration
folder_names=data.frame(folder_names=paste0("output_run_",seq(1,6)))

# Import glaicer extents
GLC_2006=raster("C:/Users/Asus/Documents/StudyProject2023S/WaSiM_setuptest/Init/glc_2006_larger_domain.asc")
# Set all NA values to 0
#values(GLC_2006)[values(GLC_2006)==NA]=0
GLC_2006[is.na(GLC_2006[])]=0

# Process the simulated snow cover maps by adding the glacier extents from 2006
# Path to the sensitivity runs 
#setwd("C:/Users/Asus/Documents/StudyProject2023S/WaSiM_setuptest")

for (n in 1:6) {
  print(n)
  # Glaicer cells have to be larger than 5 --> otherwise they are reclassified as no snow
  GLC_2006[GLC_2006==1]=10
  ssto_list = list.files(paste0("./",folder_names$folder_names[n]),full.names=TRUE,"sstodem_25_2006_large_")
  # Data frame with the mean daily fractual snow cover
  fsc_values=data.frame(date=ymd(str_extract(ssto_list,'[0-9]{4}.[0-9]{2}')),fsc_mean=0,fsc_sd=0)

  for (x in 1:length(ssto_list)) {
    print(x)
    ssto_raster=raster(ssto_list[x])
    # Set glacier to constant snow covered --> with glacier extents from 2006
    ssto_raster_glc=overlay(ssto_raster,GLC_2006,fun=function(x,y){(x+y)})
    reclass_ssto=matrix(c(0,5,0,
                          5,1000,1,
                          1000,50000,1),ncol=3, byrow=TRUE)
    ssto_raster_glc_final=reclassify(ssto_raster_glc,reclass_ssto)
    fsc_values[x,2]=round(mean(values(ssto_raster_glc_final),na.rm=T),2)
    fsc_values[x,3]=round(sd(values(ssto_raster_glc_final),na.rm=T),2)
    writeRaster(ssto_raster_glc_final,filename = paste0("./",folder_names$folder_names[n],"/Snow_maps_",fsc_values$date[x],".asc"),driver = "GeoTiff",overwrite=TRUE)
    
  } 
  write.table(fsc_values,paste0("./",folder_names$folder_names[n],"/Snow_maps_fsc.csv"),quote = F,row.names = F,sep = "\t")
}
 
