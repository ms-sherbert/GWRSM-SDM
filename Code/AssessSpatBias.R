#=== Assessment of spatio-temporal biases in swamp maire records ===#
# Written by: SM Herbert
# Written for R version 4.3.1
# This has been partially updated to use the terra and sf packages 
# (rather then the depreciated raster and sp packages)
# Last tested: 05 June 2025


#===== Preamble ======# 

rm(list=ls())

#Set up some shorthand objects for file paths
local.dir <- "D:/" #Change so that file path matches directory where repositories are cloned to your computer
local.files <- "D:/Repositories/Offline-files-GWRSM-SDM/" #file path where any gitignored files are stored locally

#Remember that working directory needs to be set to the local copy of the GWRSM-SDM repository, e.g.
setwd(paste0(local.dir,"Repositories/GWRSM-SDM"))

#---install occAssess package if required---#
#install.packages("devtools")
#require("devtools")
#devtools::install_github("robboyd/occAssess")

library(occAssess)
library(sf) #future replacement for "sp" library
library(readr)
library(raster)
require(sp)
library(terra) #future replacement for "raster" library
library(gridExtra)
library(tidyverse)
library(scales)
library(geodata)

#================= Read in data and manipulate into format required ===================#

SMrecords <- read.csv(paste0(local.files,"GWRmyrtles.csv")) # Note that this line won't work unless you have the files in your local directory
GWR <- read_sf("GWRboundary/GWRboundary2193.shp")

# Assign a plausible amount of error (in metres) to records where error not recorded
# The year 2000 used as cut-point because civilians able to access non-degraded GPS signals 
# globally from May 2nd 2000
SMrecords["Lat_Error"][is.na(SMrecords["Lat_Error"]) & SMrecords["Date_collect"] > 2000] <- 6
SMrecords["Long_Error"][is.na(SMrecords["Long_Error"]) & SMrecords["Date_collect"] > 2000] <- 6
SMrecords["Lat_Error"][is.na(SMrecords["Lat_Error"]) & SMrecords["Date_collect"] < 2001] <- 10
SMrecords["Long_Error"][is.na(SMrecords["Long_Error"]) & SMrecords["Date_collect"] < 2001] <- 10

# mask <- terra::rast(paste0(local.files,"LUCI-flood-risk/floodrisk2023.tif")) #Doesn't currently work with occAssess, commented for presumed future migration to terra package
mask <- raster::raster(paste0(local.files,"LUCI-flood-risk/floodrisk2023.tif")) # Note that this line won't work unless you have the files in your local directory

str(SMrecords) #check data type in each column in your dataframe

Error<-ifelse(SMrecords$Lat_Error > SMrecords$Long_Error, 
              SMrecords$Lat_Error, SMrecords$Long_Error) #Use whichever error larger

SMrecords<-cbind(SMrecords,Error) #append to data frame

#================ Assess spatial bias ====================#

#specify time periods for splitting data
decade_agg = list(1868:2009,2010:2022) 

SpatB<-assessSpatialBias(dat=SMrecords,
                         species="Taxon",
                         year="Date_collect",
                         x="Long_DD",
                         y="Lat_DD", 
                         spatialUncertainty="Error",
                         identifier="Group", 
                         periods=decade_agg,mask=mask,
                         degrade=TRUE) #Takes a while to run

SpatBplot<-SpatB$plot + 
  scale_x_discrete(limits = c("1", "2"),
                   labels = c("Pre-2010","2010-2022")) +
  theme_bw()

#=========== Assess environmental bias vs WorldClim data ===============#

#clim <- raster::getData("worldclim",var="bio",res=0.5,lon=175,lat=-41) #code previously used; now depreciated
clim <- geodata::worldclim_global(var="bio",path=tempdir(),res=0.5,lon=175,lat=-41) #this works but takes a long time

# delineate ROI
shp <- raster::shapefile("GWRboundary/GWRboundary2193.shp")
shp <- sp::spTransform(shp, raster::crs(clim))
#clim <- raster::crop(clim, raster::extent(shp)) 
#clim <- raster::mask(clim, shp)
clim <- terra::crop(clim, shp, mask = TRUE) #update of previous two lines (depreciated) for terra package

## extract climate data at coordinates of occurrence data 
envDat <- raster::extract(clim, SMrecords[, c("Long_DD", "Lat_DD")])


## extract background environmental data 
#backgroundEnvDat <- raster::sampleRandom(clim, size = 60, xy = F) #depreciated code
backgroundEnvDat <- terra::spatSample(clim, size = 60, method = "random", xy = FALSE)


#specify time periods
all_years = list(1868:2022) 

#This needs a bug fix - to do in v 1.0.1
envBias <- assessEnvBias(dat = SMrecords,
                         species = "Taxon",
                         x = "Long_DD",
                         y = "Lat_DD",
                         year = "Date_collect", 
                         spatialUncertainty = "Error",
                         identifier = "Group",
                         envDat = envDat[,2:20],
                         backgroundEnvDat = backgroundEnvDat,
                         xPC = 1,
                         yPC = 2,
                         periods = all_years) 

envBias$plot + theme_bw()
