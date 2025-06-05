#====== Assessment of spatial bias against common biases =====#
# Written by: SM Herbert
# Written for R version 4.3.1
# Last test: 5 June 2025

#--- Preamble ---#

rm(list=ls())

#Set up some shorthand objects for file paths
local.dir <- "D:/" #Change so that file path matches directory where repositories are cloned to your computer
local.files <- "D:/Repositories/Offline-files-GWRSM-SDM/" #file path where any gitignored files are stored locally

setwd(paste0(local.dir,"Repositories/GWRSM-SDM"))

library(sampbias)
library(raster)
require(sp)
library(terra) #future replacement for "raster" library
library(sf) #future replacement for "sp" library
library(gridExtra) #prob not needed
library(dplyr)
library(ggplot2)
library(scales) #not needed?

#--- Read in data and manipulate into format required ---#

SMrecords <- read.csv(paste0(local.files,"SMobs269.csv")) # Note that this line won't work unless you have the files in your local directory
GWR <- read_sf("GWRboundary/GWRboundary2193.shp")

#--- Bias assessment ---#

test<- sampbias::calculate_bias(x = SMrecords,
                      res = 0.01,
                      buffer = 0.5,
                      restrict_sample = NULL,
                      verbose = FALSE)

summary(test)
plot(test)

proj <- project_bias(test)
map_bias(proj)