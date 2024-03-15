#=== Assessment of spatio-temporal biases in swamp maire records ===#
# Written by: SM Herbert
# Written for R version 4.3.1


#===== Preamble ======# 

rm(list=ls())
setwd("D:/Repositories/GWRSM-SDM") #modify as needed for your computer

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

#--- Read in data and manipulate into format required ---#

SMrecords <- read.csv("SM_obs/GWRmyrtles.csv") 
GWR <- read_sf("GWRboundary/GWRboundary2193.shp")

# Assign a plausible amount of error (in metres) to records where error not recorded
# The year 2000 used as cut-point because civilians able to access non-degraded GPS signals 
# globally from May 2nd 2000
SMrecords["Lat_Error"][is.na(SMrecords["Lat_Error"]) & SMrecords["Date_collect"] > 2000] <- 6
SMrecords["Long_Error"][is.na(SMrecords["Long_Error"]) & SMrecords["Date_collect"] > 2000] <- 6
SMrecords["Lat_Error"][is.na(SMrecords["Lat_Error"]) & SMrecords["Date_collect"] < 2001] <- 10
SMrecords["Long_Error"][is.na(SMrecords["Long_Error"]) & SMrecords["Date_collect"] < 2001] <- 10

#mask <- terra::rast("GISinputs-repositories/LUCI-flood-risk/floodrisk2023.tif") 
#Doesn't currently work with occAssess, commented for presumed future migration to terra package
mask <- raster::raster("GISinputs-repositories/LUCI-flood-risk/floodrisk2023.tif")

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
                         degrade=TRUE) #NB takes a while to run

SpatBplot<-SpatB$plot + 
  scale_x_discrete(limits = c("1", "2"),
                   labels = c("Pre-2010","2010-2022")) +
  theme_bw()

#=========== Assess environmental bias vs WorldClim data ===============#

clim <- raster::getData("worldclim",var="bio",res=0.5,lon=175,lat=-41)

# delineate ROI
shp <- raster::shapefile("GWRboundary/GWRboundary2193.shp")
shp <- sp::spTransform(shp, raster::crs(clim))
clim <- raster::crop(clim, raster::extent(shp))
clim <- raster::mask(clim, shp)

## extract climate data at coordinates of occurrence data 
envDat <- raster::extract(clim, SMrecords[, c("Long_DD", "Lat_DD")])

## extract background environmental data 
backgroundEnvDat <- raster::sampleRandom(clim, size = 60,
                                         xy = F)

#specify time periods
all_years = list(1868:2022) 

envBias <- assessEnvBias(dat = SMrecords,
                         species = "Taxon",
                         x = "Long_DD",
                         y = "Lat_DD",
                         year = "Date_collect", 
                         spatialUncertainty = "Error",
                         identifier = "Group",
                         envDat = envDat,
                         backgroundEnvDat = backgroundEnvDat,
                         xPC = 1,
                         yPC = 2,
                         periods = all_years) 

envBias$plot + theme_bw()
