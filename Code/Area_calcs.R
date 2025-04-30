# Code for calculating area of suitable habitat under different scenarios
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

#=== Read in files ===#

# Study domain
clipout <-
  st_read("LCDB5-clipping-layers/LCDB5-open-water-and-rivers.shp")
GWR <-
  st_read("GWRboundary/GWRboundary2193.shp") %>%
  st_difference(y = st_union(clipout))

#Unfiltered intensity predictions
iPPMintensity <- st_read("Model-outputs/iPPM_predictions_lattice100_2025-01-22.shp")
names(iPPMintensity) #"scld_mn"

#Percentage cover waterlogged soils per 100 x 100 m pixel
prop_suitable <- rast(paste0(local.files,"Proportion_suitable_soils.tif"))

#Scenario 1 infection risk (is the binary raster for Scenario 1 - IR values that excludes 100% of myrtle rust observations)
IRS1 <- rast(paste0(local.files,"Scenario1.tif"))
IRS1 <- crop(IRS1, GWR, mask = TRUE)

#Scenario 2 infection risk (is the binary raster for Scenario 6 - IR values that excludes 75% of myrtle rust observations)
IRS2 <- rast(paste0(local.files,"Scenario6.tif"))
IRS2 <- crop(IRS2, GWR, mask = TRUE)

#Distance to road raster
EUCroad <- rast(paste0(local.files,"Euclid_distance_road/dis2road7Euc500r.tif"))%>%
  project("EPSG:2193")
EUCroad <- crop(EUCroad, GWR, mask = TRUE)
Road <- terra::ifel(EUCroad == "2663",1,0) #change values to 0 and 1

#=== Extract raster values at predicted intensity points ===#

Suitability <- terra::extract(prop_suitable,iPPMintensity)
Scenario1_IR <- terra::extract(IRS1,iPPMintensity)
Scenario2_IR <- terra::extract(IRS2,iPPMintensity)
road_pts <- terra::extract(Road,iPPMintensity)

#=== Calculations of areas under different management scenarios ===#

All <- Suitability$floodrisk2023 #Zero is not suitable
IR1 <- Suitability$floodrisk2023 * Scenario1_IR$lyr.1
IR2 <- Suitability$floodrisk2023 * Scenario2_IR$lyr.1
Accessible <- Suitability$floodrisk2023 * road_pts$Area
Accessible_IR1 <- Suitability$floodrisk2023 * road_pts$Area * Scenario1_IR$lyr.1
Accessible_IR2 <- Suitability$floodrisk2023 * road_pts$Area * Scenario2_IR$lyr.1

#Group intensity values into "high" and "low"

Intensity <- iPPMintensity$scld_mn
Intensity[Intensity < 0] <- 0 #Note that there are no zero values in this vector
Intensity[Intensity > 0] <- 1

Mgmt_scenarios <- as.data.frame(cbind(Intensity,All, IR1, IR2, Accessible, Accessible_IR1, Accessible_IR2))

Mgmt_scenarios["Intensity"][Mgmt_scenarios["Intensity"] == 0] <- "Low"
Mgmt_scenarios["Intensity"][Mgmt_scenarios["Intensity"] == 1] <- "High"

Scenarios <- pivot_longer(Mgmt_scenarios,cols=c(2:7),names_to="Scenario",values_to="Area",values_drop_na = TRUE)

Area_km2 <- Scenarios %>%
          group_by(Scenario, Intensity) %>%
          summarise("Waterlogged_area_km2" = sum(Area,na.rm=TRUE)/100,
                    "Total_area_km2" = sum(ifelse(Area > 0,1,0)/100)) 

write.csv(Area_km2,paste0("Model-outputs/Area_suitable-",run.date,".csv"))
