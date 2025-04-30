# Bootstrapping code for evaluating swamp maire presence points against LUCI model predictions
# Written by: S Herbert
# Last modified: 21/02/2024
# Written for R version 4.3.1 (S Herbert)

# Description of what this code does:

# We combined the predictions from the PPM model with the highest prediction accuracy (i.e. lowest CRPS)
# with LUCI classifications of flood accumulation following the methods of Fournier et al. (2017). 
# In brief, we compared the observed number of swamp maire observations in areas with differing LUCI 
# classifications of flood accumulation to that expected by chance from a random distribution.

# Step 1: Bootstrap 268 points 10,000 times and count the number of points sampled in each LUCI category.
# To obtain expected numbers of swamp maire records in each LUCI flood accumulation category.

# Step 2: Classify habitats as suitable (coded 1) if the number of observed points 
# is in the 97.5th percentile of the bootstrapped distribution. 
# Classify habitats as unsuitable (coded -1) if the observed number of points is in the 
# 2.5th percentile of the bootstrapped distribution. 
# Classify habitats neutral if the number of observed points fell between the 2.5th and 97.5th percentiles.  

# Step 3: Construct a binary spatial filter from the LUCI layer (1 = suitable and neutral, 0 = unsuitable) 
# aggregate (to a value between 0 and 1 by mean suitability) 100 x 100 m PPM resolution to represent 
# the area of suitable habitat on a scale between 
# multiply the iPPM predictions by this to create a weighted suitability score


#=== Preamble - use bits as needed for your computer ===#

rm(list=ls())

local.dir <- "D:/"
setwd(paste0(local.dir,"Repositories/GWRSM-SDM")) #just to make sure your working directory is the cloned repository

local.files <- "D:/Repositories/Offline-files-GWRSM-SDM/" #file path where any gitignored files are stored locally
run.date <- as.character(Sys.Date())

#Load required packages
#I think the inlabru / INLA / PointedSDMs packages will be need to recognise the BruSDM outputs
library(PointedSDMs) #version 1.3
require(inlabru)     #version 2.10.1
require(fmesher)
require(R6)
library(sf)
library(INLA)        #version 23.04.24 (HR Lai's computer) or version 23.09.09 (VUW PC)
library(sf)
library(terra)
library(tidyverse)
library(tidyterra)

#=== Read in Data sources ===#

# Study domain
clipout <-
  st_read("LCDB5-clipping-layers/LCDB5-open-water-and-rivers.shp")
GWR <-
  st_read("GWRboundary/GWRboundary2193.shp") %>%
  st_difference(y = st_union(clipout))

#iPPM predictions

iPPMintensity <- st_read("Model-outputs/iPPM_predictions_lattice100_2025-01-22.shp")

#LUCI flood classifications
flood <- rast(paste0(local.files,"LUCI-flood-risk/floodrisk2023.tif"))

#Observed SYZmai presence-only records
#(pooled from NVS, iNaturalist, extra observation data, and herbaria)
obs <-
  read_csv(paste0(local.files,"SM-PO-PA-GWR-full.csv")) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = "+proj=longlat +ellips=WGS84") %>%
  st_transform(crs = st_crs(GWR)) %>%
  st_intersection(y = GWR)

obsSpatVec <- vect(as(obs, "Spatial"))

#=== Step 1: Bootstrap 268 points inside GWR expected from a random distribution, 10,000 times ===#

Nruns <- 10000
datalist = vector(mode = "list", length = Nruns)

for(i in 1:Nruns){
  #Bootstrap 268 points
  rand_points<- st_sample(GWR,
                          size = 268,
                          type = "random",
                          exact = TRUE,
                          warn_if_not_integer = TRUE,
                          by_polygon = FALSE,
                          progress = TRUE,
                          force = FALSE
                          )
  #Extract LUCI category at point positions
  RPts <- vect(as(rand_points, "Spatial"))
  LUCIcount_exp <- terra::extract(flood, RPts)  
  #Count how many points per 
  Count_expected <- LUCIcount_exp %>% count(floodrisk2023)
  Run_no<-rep(i,times=length(Count_expected$floodrisk2023))
  Count_expected<-cbind(Count_expected,Run_no)
  datalist[[i]] <- Count_expected  
  }

#If we want to see a plot of random points (it'll plot the points from the final bootstrap) vs observed points
plot(rand_points,pch=16,col="blue")
points(obsSpatVec,pch=16)
lines(GWR)

#=== Step 2: Habitat classification as suitable / unsuitable ===#

# Step 2: Classify habitats as suitable (coded 1) if the number of observed points 
# is in the 97.5th percentile of the bootstrapped distribution. 
# Classify habitats as unsuitable (coded -1) if the observed number of points is in the 
# 2.5th percentile of the bootstrapped distribution. 
# Classify habitats neutral if the number of observed points fell between the 2.5th and 97.5th percentiles.  

#Calculate percentiles of expected counts

big_data = do.call(rbind, datalist)

ExpQuantiles <- big_data %>% 
                group_by(floodrisk2023) %>%  
                  summarise(quantile = c("2.5%","97.5%"),
                  n = quantile(n, c(0.025, 0.975))
                  )

#Calculate observed counts
LUCIcount_obs <- terra::extract(flood, obsSpatVec)  
Count_obs <- LUCIcount_obs %>% count(floodrisk2023)

#If we want to see a histogram of (our last bootstrap of) expected vs observed counts

Observed <-  ggplot(data = Count_obs,aes(x=floodrisk2023,y=n)) +
  geom_bar(stat="identity") +
  ylab("Number of points") +
  xlab("Soil waterlogging category") +
  ggtitle("Observed") +
  theme_bw()

Expected <-  ggplot(data = Count_expected,aes(x=floodrisk2023,y=n)) +
  geom_bar(stat="identity") +
  ylab("Number of points") +
  xlab("Soil waterlogging category") +
  ggtitle("Expected") +
  theme_bw()

multiplot(Observed, Expected)

#View Count_obs and ExpQuantiles tibbles and just evaluate suitability manually
Count_obs
ExpQuantiles

# LUCI 1: Observed points in < 2.5% quantile - unsuitable - code 0
# LUCI 2: Observed points in < 2.5% quantile - unsuitable - code 0
# LUCI 3: Observed points in < 2.5% quantile - unsuitable - code 0
# LUCI 4: Observed points in > 97.5% quantile - suitable - code 1
# LUCI 5: Observed points in > 97.5% quantile - suitable - code 1


#=== Step 3: Assign suitability scores and aggregate to a 100 x 100m resolution ===#

FloodSuit <- flood

# Replace flood risk 2023 values with binary suitability index
FloodSuit <- subst(FloodSuit, 1, 0)
FloodSuit <- subst(FloodSuit, 2, 0)
FloodSuit <- subst(FloodSuit, 3, 0)
FloodSuit <- subst(FloodSuit, 4, 1)
FloodSuit <- subst(FloodSuit, 5, 1)

FloodSuit$floodrisk2023 #check it worked

FloodSuit_100 <- aggregate(FloodSuit, fact=20)
res(FloodSuit_100) #Double-check that the resolution is correct
plot(FloodSuit_100)

FloodSuit_100_cropped <- crop(FloodSuit_100, GWR, mask = TRUE)
writeRaster(FloodSuit_100_cropped,"D:/Repositories/GWRSM-SDM/Proportion_suitable_soils.tif", overwrite =TRUE, filetype = "GTiff")

Suitability <- terra::extract(FloodSuit_100,iPPMintensity) #interpolate predicted iPPM point locations with FloodSuit_100 raster to get multiplier values
Filtered_suitability <- iPPMintensity$scld_mn * Suitability
Filtered_suitability_mean <- iPPMintensity$mean * Suitability
Filtered_suitability_scaled <- scale(Filtered_suitability_mean$floodrisk2023)
iPPMintensity <- cbind(iPPMintensity,Filtered_suitability$floodrisk2023,Filtered_suitability_mean$floodrisk2023,Filtered_suitability_scaled)

Filtered_points <- ggplot() +
  geom_sf(data=iPPMintensity, aes(geometry = geometry, color = Filtered_suitability.floodrisk2023)) +
  geom_sf(data = GWR, alpha = 0.1) +
  #geom_sf(data = subset(PAobs, NPres == 0), color = "black") +
  #geom_sf(data = POobs, color = "grey70",pch=15) +
  #geom_sf(data = subset(PAobs, NPres == 1), color = "grey70") +
  ggtitle("Filtered predicted mean intensity (scaled)") +
  scale_color_viridis_c(name = "Filtered predicted mean intensity (scaled)") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

#Export projected values to a shapefile for analysis in ArcGIS or QGIS
st_write(iPPMintensity,paste0("Model-outputs/filtered_predictions_100_",run.date,".shp"),overwrite=TRUE)

#Export map of filtered projections

ggsave(Filtered_points,
       filename=paste0("Model-outputs/Filteredpreds100m_",run.date,".png"),
       device="png",
       height=10, width = 16, units = "cm",dpi="print")

#Calculate how many scaled points were predicted to be greater than the average 

length(which(iPPMintensity$scld_mn > 0))
