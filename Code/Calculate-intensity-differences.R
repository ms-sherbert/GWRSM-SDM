# Create delta figure for supplementary material
# Written for R version 4.3.1

#=== Preamble - use bits as needed for your computer ===#

rm(list=ls())

local.dir <- "D:/"
setwd(paste0(local.dir,"Repositories/GWRSM-SDM")) #just to make sure your working directory is the cloned repository

local.files <- "D:/Repositories/Offline-files-GWRSM-SDM/" #file path where any gitignored files are stored locally
run.date <- as.character(Sys.Date())

library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(gdalraster)


#=== Read in files ===#

# Study domain
clipout <-
  st_read("LCDB5-clipping-layers/LCDB5-open-water-and-rivers.shp")
GWR <-
  st_read("GWRboundary/GWRboundary2193.shp") %>%
  st_difference(y = st_union(clipout))

#Unfiltered intensity predictions
iPPMintensity <- st_read(paste0(local.files,"iPPM_predictions_lattice100_2025-01-22.shp"))
POPPMintensity <- st_read(paste0(local.files,"PO-PPM_predictions_100_2024-08-19.shp"))

Delta_mean <- POPPMintensity$mean-iPPMintensity$mean 

iPPMintensity$Delta <- Delta_mean

st_write(iPPMintensity,paste0(local.files,"Deltas.shp")) 

gdalraster::rasterize(src_dsn=paste0(local.files,"Deltas.shp"),
                      dstfile=paste0(local.files,"Deltas.tif"),
                      layer="Delta",tr=c(100,100),dstnodata=0,fmt="GTiff")


  