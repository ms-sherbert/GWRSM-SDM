#====== Assessment of spatial bias against common biases =====#
# Written by: S Herbert
# Written for R version 4.3.1


#--- Preamble ---#

rm(list=ls())
setwd("D:/Repositories/GWRSM-SDM") #or wherever the repository is cloned

library(sampbias)
library(raster)
require(sp)
library(terra) #future replacement for "raster" library
library(sf) #future replacement for "sp" library
#library(gridExtra) #prob not needed
library(dplyr)
library(ggplot2)
#library(scales) #not needed?

#--- Read in data and manipulate into format required ---#

SMrecords <- read.csv("SM_obs_GBIF/SMobs269.csv") 
GWR <- read_sf("GWRboundary/GWRboundary2193.shp")

#---

test<- sampbias::calculate_bias(x = SMrecords,
                      res = 0.01,
                      buffer = 0.5,
                      restrict_sample = NULL,
                      verbose = FALSE)

summary(test)
plot(test)

proj <- project_bias(test)
map_bias(proj)