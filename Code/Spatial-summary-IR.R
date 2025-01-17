# Read in infection risk rasters downloaded from the NIWA website 
# and use to summarise spatial variation in risk between August 2020 and January 2025
# Code originally written by Stephanie Tomscha, June 2022
# Process updated to use terra and sf packages by Sarah Herbert, January 2025
# Also changed to spatiotemporal summary scenarios of infection risk
# Written for R version 4.3.1

#=== Preamble - use bits as needed for your computer ===#

rm(list=ls())

library(terra)
library(sf)
library(tidyverse)
library(tidyterra)

#=== Required functions and definitions ===#

counter <- function(x) length(which(x > 0.5)) / length(which(x>=0)) 
#calculates the proportion of weeks in which mean IR > 0.5

NZTM<-"+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"

#=== Read in data required and format ===#

current.list <- list.files(path="Z:/BioProtectAotearoa/New_meanIR_download", 
                           pattern =".asc$", full.names=TRUE)
s <- terra::rast(current.list) #create raster stack from all asci files
crs(s) <- "epsg:2193"

GWR <- st_read("D:/Repositories/GWRSM-SDM/GWRboundary/GWRboundary2193.shp")

scrop <- crop(s,GWR, mask=TRUE) #crop raster stack to GWR polygon

MRobs <- read_csv("D:/Repositories/GWRSM-SDM/MRobs/iNaturalist-MRobs-Aug2020-Jan2025.csv") %>%
         st_as_sf(coords=c("longitude","latitude"),
                  crs = "+proj=longlat +ellips=WGS84") 

MRobsV <- vect(MRobs) %>%
          project("EPSG:2193")

#=== Summarise weekly mean infection risk over c. 4 year period ===# 

#smax<-max(scrop)
#smean<-mean(scrop)
scount<-app(scrop,counter) #Preferred method - proportion of weeks where predicted IR > 50%
  
#Optional - code to plot summary outputs to check
#plot(scount)
#polys(GWR, border="black", lwd=1, lty=1, alpha=1) #Show GWR boundary
#points(MRobsV, col= "black", pch=16) #Add myrtle rust observations

#=== Extract predicted IR thresholds based on percentiles of observed MR cases ===#

IRvMRobs <- terra::extract(scount, MRobsV, method="simple")

Scenarios <- quantile(IRvMRobs$lyr.1, probs = c(0, 0.05, 0.1, 0.15),na.rm=TRUE)

Scenario1 <- terra::ifel(scount < 0.01716738,1,0) #acceptable proportion of IR50-weeks at 0th percentile of MR obs
Scenario2 <- terra::ifel(scount < 0.02145923,1,0) #acceptable proportion of IR50-weeks at 5th percentile of MR obs
Scenario3 <- terra::ifel(scount < 0.02769671,1,0) #acceptable proportion of IR50-weeks at 10th percentile of MR obs
Scenario4 <- terra::ifel(scount < 0.03004292,1,0) #acceptable proportion of IR50-weeks at 15th percentile of MR obs

#plot scenarios against myrtle rust observation points
par(mfrow=c(2,2))

plot(Scenario1,main="Scenario 1: 0th percentile of myrtle rust observations")
polys(GWR, border="black", lwd=1, lty=1, alpha=1)
points(MRobsV, col= "black", pch=16)

plot(Scenario2,main="Scenario 2: 5th percentile of myrtle rust observations")
polys(GWR, border="black", lwd=1, lty=1, alpha=1)
points(MRobsV, col= "black", pch=16)

plot(Scenario3,main="Scenario 3: 10th percentile of myrtle rust observations")
polys(GWR, border="black", lwd=1, lty=1, alpha=1)
points(MRobsV, col= "black", pch=16)

plot(Scenario4,main="Scenario 4: 15th percentile of myrtle rust observations")
polys(GWR, border="black", lwd=1, lty=1, alpha=1)
points(MRobsV, col= "black", pch=16)


#Export scount and binarised IR scenario rasters

writeRaster(scount,"D:/Repositories/GWRSM-SDM/scount.tif", overwrite =TRUE, filetype = "GTiff")
writeRaster(Scenario1,"D:/Repositories/GWRSM-SDM/Scenario1.tif", overwrite =TRUE, filetype = "GTiff")
writeRaster(Scenario2,"D:/Repositories/GWRSM-SDM/Scenario2.tif", overwrite =TRUE, filetype = "GTiff")
writeRaster(Scenario3,"D:/Repositories/GWRSM-SDM/Scenario3.tif", overwrite =TRUE, filetype = "GTiff")
writeRaster(Scenario4,"D:/Repositories/GWRSM-SDM/Scenario4.tif", overwrite =TRUE, filetype = "GTiff")