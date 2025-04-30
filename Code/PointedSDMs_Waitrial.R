# Modification of Hao Ran's inlabru script to explore the feasibility of integrated SDM using point process model
# Written by: HR Lai and SM Herbert
# Written for R version 4.3.1 (SMH) and 4.2.2 (HRL)


#=== Preamble - use bits as needed for your computer ===#

rm(list=ls())

#Set up some shorthand objects for file paths
local.dir <- "D:/" #Change so that file path matches directory where repositories are cloned to your computer
local.files <- "D:/Repositories/Offline-files-GWRSM-SDM/" #file path where any gitignored files are stored locally

#Remember that working directory needs to be set to the local copy of the GWRSM-SDM repository, e.g.
setwd(paste0(local.dir,"Repositories/GWRSM-SDM"))

run.date <- as.character(Sys.Date())

#Use if needed to install PointedSDMs package

#If running R v4.3.2 or more recent; presumably install from CRAN works
#install.packages('PointedSDMs')

#If running older version of R (e.g. 4.3.1)
#devtools::install_github("PhilipMostert/PointedSDMs")

#install.packages("INLA",repos=c(getOption("repos"),
#                INLA="https://inla.r-inla-download.org/R/stable"),
#                dep=TRUE)

#We are having some issues with terra if you try to read the saved bru.sdm files from .rds
#This version is the most stable:
#packageurl <- "http://cran.r-project.org/src/contrib/Archive/terra/terra_1.7-46.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

#Load required packages
library(PointedSDMs) #version 1.3
require(inlabru)     #version 2.10.1
require(fmesher)
require(R6)
library(sf)
library(INLA)        #version 23.04.24 (HR Lai's computer) or version 23.09.09 (VUW PC)
library(terra)       #version 1.7-46
library(tidyverse)
library(scoringRules)
library(tidyterra)
source("Code/utils.R")

#=== Read in Data sources ===#

# Study domain
clipout <-
  st_read("LCDB5-clipping-layers/LCDB5-open-water-and-rivers.shp")
GWR <-
  st_read("GWRboundary/GWRboundary2193.shp") %>%
  st_difference(y = st_union(clipout))

# Covariates
wet1 <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/topo_wetness.tif"))
wet1 <- crop(wet1, GWR, mask = TRUE)
wet <- scale(wet1)

drain1 <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/soil_drainage.tif"))
drain1 <- crop(drain1, GWR, mask = TRUE)
drain <- scale(drain1)

Tmin1 <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/temp_minColdMonth.tif"))
Tmin1 <- crop(Tmin1, GWR, mask = TRUE)
Tmin <- scale(Tmin1)

precip1 <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/precip_warmQtr.tif"))
precip1 <- crop(precip1, GWR, mask = TRUE)
precip <- scale(precip1)

humid1 <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/humidity_meanAnn.tif"))
humid1 <- crop(humid1, GWR, mask = TRUE)
humid <- scale(humid1)

solar1 <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/solRad_winter.tif"))
solar1 <- crop(solar1, GWR, mask = TRUE)
solar <- scale(solar1)

Trange1 <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/temp_annRange.tif"))
Trange1 <- crop(Trange1, GWR, mask = TRUE)
Trange <- scale(Trange1)

road1 <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/distance_road.tif"))
road1 <- crop(road1, GWR, mask = TRUE)
road <- scale(road1) #type in assigned name (e.g. 'road)' to check sf details

# SYZmai records

obs_all <-
  read_csv(paste0(local.files,"SM-PO-PA-GWR-full.csv")) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = "+proj=longlat +ellips=WGS84") %>%
  st_transform(crs = st_crs(GWR)) %>%
  st_intersection(y = GWR)

PAobs<-dplyr::filter(obs_all, Source == "NVS")
POobs<-dplyr::filter(obs_all, Source != "NVS")
#Take out the column that provides presence/absence in the presence-only data set or the ISDM model gets confused
drops <- c("NPres")
POobs<- POobs[ , !(names(POobs) %in% drops)]

PAobs <- PAobs %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = "+proj=longlat +ellips=WGS84") %>%
  st_transform(crs = st_crs(GWR)) %>%
  st_intersection(y = GWR)

POobs <- POobs %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = "+proj=longlat +ellips=WGS84") %>%
  st_transform(crs = st_crs(GWR)) %>%
  st_intersection(y = GWR)

#Management scenarios

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

#=== Graph data sources for Fig S4 ===#

GWRplot <- as(GWR, "Spatial")
GWRplot <- vect(GWRplot)

par(mfrow=c(4,2))
terra::plot(Tmin1,xlab="Easting",ylab="Northing",main = "Min temperature coldest month (deg C)")
terra::polys(GWRplot, border="black", lwd=1, lty=1, alpha=1)
terra::plot(precip1,xlab="Easting",ylab="Northing",main = "Precipitation warmest quarter (mm)")
terra::polys(GWRplot, border="black", lwd=1, lty=1, alpha=1)
terra::plot(humid1,xlab="Easting",ylab="Northing",main = "Mean annual humidity (%)")
terra::polys(GWRplot, border="black", lwd=1, lty=1, alpha=1)
terra::plot(solar1,xlab="Easting",ylab="Northing",main = "Solar radiation in winter (MJ per sq.m per day)")
terra::polys(GWRplot, border="black", lwd=1, lty=1, alpha=1)
terra::plot(Trange1,xlab="Easting",ylab="Northing",main = "Annual temperature range (deg C)")
terra::polys(GWRplot, border="black", lwd=1, lty=1, alpha=1)
terra::plot(drain1,xlab="Easting",ylab="Northing",main = "Soil drainage capacity")
terra::polys(GWRplot, border="black", lwd=1, lty=1, alpha=1)
terra::plot(wet1,xlab="Easting",ylab="Northing",main = "Topographic wetness index")
terra::polys(GWRplot, border="black", lwd=1, lty=1, alpha=1)
terra::plot(road1,xlab="Easting",ylab="Northing",main = "Distance from nearest road (m)")
terra::polys(GWRplot, border="black", lwd=1, lty=1, alpha=1)
par(mfrow=c(1,1)) #Re-set partitions in plotting space to a single panel

#=== Model preparation ===#

# Mesh
# following https://rpubs.com/jafet089/886687
# also https://haakonbakkagit.github.io/btopic104.html

coords <- st_coordinates(obs_all)
domain <- as(GWR, "Spatial") %>% inla.sp2segment()

max.range <- max(apply(coords, 2, function(x) diff(range(x))))
max.edge <- max.range / (3*5)
bound.outer <-  max.range / 3
mesh1 <- inla.mesh.2d(
  boundary = domain,
  loc = coords,
  max.edge = c(1, 5) * max.edge,
  offset = c(max.edge, bound.outer),
  cutoff = 500,   
  crs = st_crs(GWR)
)

mesh1$n
# plot(mesh1)
ggplot() +
  gg(data = GWR, fill = "darkgreen", colour = NA, alpha = 0.2) +
  gg(mesh1) +
  gg(data = subset(PAobs, NPres== 1), pch = 21, fill = "white") +
  gg(data = subset(PAobs, NPres == 0), pch = 21, fill = "black") +
  gg(data = POobs, pch = 21, fill = "purple") +
  ggtitle(label = "B") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

proj <- sp::CRS("+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000
                  +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m
                +no_defs +type=crs")

# Stack fixed factor covariates together into one raster object
# Two model options are provided here, but we only use the global model formulation for predictions

covariates <- c(wet,precip,solar,Tmin,Trange,humid,road,drain) #global model
#covariates_minimal <-c(solar,drain,road) #minimal model containing only significant effects in global model


#================= Modelling ==================#

# Refer to: https://cran.r-project.org/web/packages/PointedSDMs/vignettes/Solitary_tinamou.html
# And Mostert and O'Hara (2023)

spatial_data <- intModel(POobs,PAobs,
                         Coordinates =c("X","Y"),
                         Projection = proj,
                         Mesh = mesh1,
                         responsePA = 'NPres',
                         Boundary = GWR,
                         spatialCovariates = covariates)

#If we want to specify priors
spatial_data$priorsFixed("topo_wetness", mean.linear = 0, prec.linear = 1/2)
spatial_data$priorsFixed("precip_warmQtr", mean.linear = 0, prec.linear = 1/2)
spatial_data$priorsFixed("solRad_winter", mean.linear = 0, prec.linear = 1/2)
spatial_data$priorsFixed("temp_minColdMonth", mean.linear = 0, prec.linear = 1/2)
spatial_data$priorsFixed("temp_annRange", mean.linear = 0, prec.linear = 1/2)
spatial_data$priorsFixed("humidity_meanAnn", mean.linear = 0, prec.linear = 1/2)
spatial_data$priorsFixed("distance_road", mean.linear = 0, prec.linear = 1/2)
spatial_data$priorsFixed("soil_drainage", mean.linear = 0, prec.linear = 1/2)
spatial_data$specifySpatial(sharedSpatial = TRUE, prior.sigma = c(1, 0.5), prior.range = c(1000, 0.25))

dataset_plot <- spatial_data$plot(Boundary=FALSE) +
                geom_sf(data = GWR, alpha = 0.1) +
                theme_bw()

spat_model <- fitISDM(data = spatial_data,
                      options = list(control.inla = list(int.strategy='eb'))
                      )

spat_model

#Export model object as an .rds file so we can make more predictions later if needed
saveRDS(spat_model,paste0("Model-outputs/iPPM_",run.date,".rds"))
#spat_model <- read_rds("Model-outputs/iPPM_2024-08-12.rds") #e.g. code if need to re-load

#Export model coefficients for fixed effects
write.csv(spat_model$summary.fixed,paste0("Model-outputs/iPPM_fixed_effects_",run.date,".csv"))

#=== Calculate prediction score ===#

# Integrated predictions
# First we partition the landscape into grids of some resolution

B1 <- partition(samplers = GWR, resolution = 50*50)
# note that resolution defines the area of the grid squares in square metres
# plot(B1) - this resolution is 2500 square metres
# here's some code to export for checking resolution if needed
# st_write(B1,"GISinputs-repositories/B1_sampler.shp",append=FALSE)

# add the total observed number of occurrence points
B1$NPres <- lengths(st_intersects(B1, filter(obs_all, NPres != 0)))

# predict the total number of occurrences (counts) per grid
# As per our prediction of total abundance in the GWR,
# We leave the following line out of this chunk:
# PAobs_intercept +
# This is so that the model only predicts the number of points where the species is present

Lambda <- predict(
  spat_model,
  fm_int(domain = mesh1, sampler = B1),  # this is to integrate the prediction over grid area
  formula = ~ tapply(weight * exp(shared_spatial +
                                    POobs_intercept +
                                    topo_wetness +
                                    precip_warmQtr +
                                    solRad_winter +
                                    temp_minColdMonth +
                                    temp_annRange +
                                    humidity_meanAnn +
                                    distance_road + #For CRPS we want road back in
                                    soil_drainage),
                     .block,
                     sum)
)

abun_total <-
  Lambda$predictions %>%
  rownames_to_column(".block")

# add predicted total counts to the grid object
B1 <-
  B1 %>%
  rownames_to_column(".block") %>%
  left_join(abun_total)

# some plots
plot(select(B1, "mean"))
plot(select(B1, "median"))
with(B1, plot(NPres, median)); abline(0, 1)

# calculate CRPS
# see also https://inlabru-org.github.io/inlabru/articles/prediction_scores.html
CRPS_block <- crps_pois(B1$NPres, B1$median)
hist(CRPS_block, breaks = 50)   # we want these values to be a close to zero as possible

median(CRPS_block, na.rm = TRUE)
with(B1, plot(median, CRPS_block))

#==== Make predictions ====#
# The fixed factor covariate 'distance from road' is conditioned out of the predictive model (Warton et al. 2013)
# This is to remove the directional bias from our predictions
# i.e. make predictions over the whole range that would be true if every point in our space is equidistant from a road

projections <- predict(spat_model,
                       data=fm_pixels(mesh=mesh1, mask = GWR, dims = c(1445,1021)), #dims = c(722.5,510.5)) is for 200x200 projection
                       formula = ~
                         shared_spatial +
                         POobs_intercept +
                         # PAobs_intercept +  # I think we need to mute this to only use the Presence-only intercept
                         topo_wetness +
                         precip_warmQtr +
                         solRad_winter +
                         temp_minColdMonth +
                         temp_annRange +
                         humidity_meanAnn +
                         soil_drainage,
                       # predictor = TRUE,
                       # spatial = TRUE,
                       fun='linear',
                       n.samples = 1000)

projection_means <- plot(projections, plot=FALSE) #create ggplot object from projections
#Or to make a basic plot, just use plot(projections, plot = TRUE)

#=== Plot predictions for Figure 2 ===#
# Code for a few options for this figure are provided
# Previous formula used for converting predicted local abundance to predicted site occupancy
# projected_occupancy <- 1-exp(-D) #Relationship of occupancy and abundance (parameter D) after Caughley (1977)

#Plot predicted mean point intensity on linear predictor scale
LinearPreds <- projection_means$predictions$mean +
                geom_sf(data = GWR, alpha = 0.1) +
                #geom_sf(data = subset(PAobs, NPres == 0), color = "black",show.legend=TRUE) +
                #geom_sf(data = POobs, color = "grey70",pch=15,show.legend="point") +
                #geom_sf(data = subset(PAobs, NPres == 1), color = "grey70") +
                ggtitle("B") +
                scale_color_viridis_c(name = "Predicted mean intensity") +
                xlab("Longitude") +
                ylab("Latitude") +
                theme_bw()

#Plot scaled predicted mean point intensity on linear predictor scale
#see Morera-Pujol et al. (2023) for rationale behind scaling

scaled_mean <- scale(projections$predictions$mean)
scaled_df<-cbind(scaled_mean,projections$predictions)

ScaledPreds <- ggplot() +
                geom_sf(data=scaled_df, aes(geometry = geometry, color = scaled_mean)) +
                geom_sf(data = GWR, alpha = 0.1) +
                #geom_sf(data = subset(PAobs, NPres == 0), color = "black") +
                #geom_sf(data = POobs, color = "grey70",pch=15) +
                #geom_sf(data = subset(PAobs, NPres == 1), color = "grey70") +
                ggtitle("C") +
                scale_color_viridis_c(name = "Predicted mean intensity (scaled)") +
                xlab("Longitude") +
                ylab("Latitude") +
                theme_bw()

#Plot standard deviation of mean point intensity (linear predictor scale)
SdPreds <- ggplot() +
            geom_sf(data=scaled_df, aes(geometry = geometry, color = sd)) +
            geom_sf(data = GWR, alpha = 0.1) +
            #geom_sf(data = subset(PAobs, NPres == 0), color = "black") +
            #geom_sf(data = POobs, color = "grey70",pch=15) +
            #geom_sf(data = subset(PAobs, NPres == 1), color = "grey70") +
            ggtitle("Standard deviation of predicted mean intensity - integrated PPM") +
            scale_color_viridis_c(name = "Standard deviation") +
            xlab("Longitude") +
            ylab("Latitude") +
            theme_bw()

#Export projected values to a csv file
write.csv(scaled_df,paste0("Model-outputs/iPPM_scaled_predictions100_",run.date,".csv"))
#Export projected values to a shapefile for analysis in ArcGIS or QGIS
st_write(scaled_df,paste0("Model-outputs/iPPM_predictions_lattice100_",run.date,".shp"), append=FALSE)

#Export maps of projections

ggsave(LinearPreds,
       filename=paste0("Model-outputs/LinearPreds100m",run.date,".png"),
       device="png",
       height=10, width = 16, units = "cm",dpi="print")

ggsave(ScaledPreds,
       filename=paste0("Model-outputs/ScaledPreds100m",run.date,".png"),
       device="png",
       height=10, width = 16, units = "cm",dpi="print")

ggsave(SdPreds,
       filename=paste0("Model-outputs/SdPreds100m",run.date,".png"),
       device="png",
       height=10, width = 16, units = "cm",dpi="print")

#Calculate how many scaled points were predicted to be greater than the average 
length(which(scaled_df$scaled_mean > 0))


#=== Calculate total swamp maire abundance in the GWR ===#
# See also https://inlabru-org.github.io/inlabru/articles/2d_lgcp_sf.html#estimating-abundance
# We leave the following line out of this code chunk:
# PAobs_intercept +
# This is so that the model only predicts the number of points where the species is present

Lambda_total <- predict(
  spat_model,
  fm_int(mesh1, GWR),
  formula = ~ sum(weight * exp(
    shared_spatial +
    POobs_intercept +
    (distance_road * 0 + minmax(road)[1]) +
    topo_wetness +
    precip_warmQtr +
    solRad_winter +
    temp_minColdMonth +
    temp_annRange +
    humidity_meanAnn +
    soil_drainage)),
  fun='linear',
  n.samples = 1000
)
Lambda_total

#write.csv(Lambda_total$predictions,paste0("Model-outputs/iPPM_totabundance_predictions",run.date,".csv")

# Generate posterior abundance distribution
# The 95% CrI of Lambda_total is used to define a broader range of plugin values
# to constrain posterior values of N
# We could also use a value just lower than the number of 268 observations to constrain the lower tail
# Because we know we won't have fewer than this number of trees
Trees <- predict(
  spat_model,
  fm_int(mesh1, GWR),
  formula = ~ data.frame(
    N = seq(110, 3300, 10), #original run used seq(100,1000,10)
    dpois(seq(110, 3300, 10),
          lambda = sum(weight * exp(
              shared_spatial +
              POobs_intercept +
              topo_wetness +
              precip_warmQtr +
              solRad_winter +
              temp_minColdMonth +
              temp_annRange +
              humidity_meanAnn +
              soil_drainage))
    )
  ),
  fun='linear',
  n.samples = 4000
)

#Get quantiles of the posterior abundance distribution
inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x = Trees$predictions$N, y = Trees$predictions$mean))

#Get mean of the posterior abundance distribution
inla.emarginal(identity, marginal = list(x = Trees$predictions$N, y = Trees$predictions$mean))

#Plot posterior and plugin abundance distributions for comparison
Trees$predictions$plugin_estimate <- dpois(Trees$predictions$N, lambda = Lambda_total$predictions$median)

PostN <- ggplot(data = Trees$predictions) +
          geom_line(aes(x = N, y = mean, colour = "Posterior")) +
          geom_line(aes(x = N, y = plugin_estimate, colour = "Plugin")) +
          theme_bw()

#=== Determine land area suitable for protection ===#

# Partition the landscape into grids of 1 x 1 km resolution to match the MRPM
# note that resolution defines the area of the grid squares in square metres
B2 <- partition(samplers = GWR, resolution = 1000)
#plot(B2)
#st_write(B2,"GISinputs-repositories/B2_sampler.shp",append=FALSE)

# add the total observed number of occurrence points
B2$NPres <- lengths(st_intersects(B2, filter(obs_all, NPres != 0)))

Lambda <- predict(
  spat_model,
  fm_int(domain = mesh1, sampler = B2),  # this is to integrate the prediction over 1 x 1 km grid area
  formula = ~ tapply(weight * exp(shared_spatial +
                                    POobs_intercept +
                                    topo_wetness +
                                    precip_warmQtr +
                                    solRad_winter +
                                    temp_minColdMonth +
                                    temp_annRange +
                                    humidity_meanAnn +
                                    (distance_road * 0 + minmax(road)[1]) +
                                    soil_drainage),
                     .block,
                     sum)
)

# Match the myrtle rust risk and road accessibility values to the shapefile B2
B2SpatVec <- vect(as(B2, "Spatial"))
B2Pts <- centroids(B2SpatVec, inside=FALSE) #currently a workaround to account for imperfect alignment of the grids

#This also works to extract the mean MRPM IR value within each grid cell:
#B2IR <- terra::extract(IR, B2SpatVec,fun=mean)

B2IRS1 <- terra::extract(IRS1, B2Pts) #values in lyr.1 column are already 0 = unacceptably high risk, 1 = acceptable risk
#Change values to binary multipliers (0 = unacceptably high risk, 1 = acceptable risk)
#B2IR["lyr.1"][B2IR["lyr.1"] == 100] <- 1
#B2IR["lyr.1"][B2IR["lyr.1"] == 128] <- NA

B2IRS2 <- terra::extract(IRS2, B2Pts) #values in lyr.1 column are already 0 = unacceptably high risk, 1 = acceptable risk 

B2Rd <- terra::extract(EUCroad, B2Pts)
Access <- as.numeric(B2Rd$Area)
B2Rd <- cbind(B2Rd,Access)
#Change values to binary multipliers (0 = inaccessible, 1 = accessible)
B2Rd["Access"][B2Rd["Access"] == 1] <- 0
B2Rd["Access"][B2Rd["Access"] == 2] <- 1

B2 <-cbind(B2,B2IRS1$lyr.1,B2IRS2$lyr.1,B2Rd$Access) #Bind scenarios to B2

# join B2 to predicted abundances
abun_total <-
  B2 %>%
  rownames_to_column(".block") %>%
  left_join(
    Lambda$predictions %>%
      rownames_to_column(".block")
      )

# Scale and write the abundance per 1 x 1 km cell to shapefile
scaled_mean_1km <- scale(abun_total$mean)
abun_total<-cbind(abun_total,scaled_mean_1km)
st_write(abun_total,"Model-outputs/Abundance1km_2.shp") 
plot(abun_total[20]) 

# construct cumulative curve in six ways
cumsum_list <- list(abun_total, abun_total, abun_total, abun_total, abun_total, abun_total)

# 1. accumulation by abundance only
cumsum_list[[1]] <-
  cumsum_list[[1]] %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "All")

# 2. accumulation by abundance in accessible areas only
cumsum_list[[2]] <-
  cumsum_list[[2]] %>%
  filter(B2Rd.Access == 1) %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "Accessible")

# 3. accumulation by abundance in Scenario 1 low myrtle rust risk areas only
cumsum_list[[3]] <-
  cumsum_list[[3]] %>%
  filter(B2IRS1.lyr.1 == 1) %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "Low infection risk Scenario 1")

# 4. accumulation by abundance in Scenario 2 low myrtle rust risk areas only
cumsum_list[[4]] <-
  cumsum_list[[4]] %>%
  filter(B2IRS2.lyr.1 == 1) %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "Low infection risk Scenario 2")

# 5. accumulation by abundance in low myrtle rust risk Scenario 1 and accessible areas
cumsum_list[[5]] <-
  cumsum_list[[5]] %>%
  filter(B2IRS1.lyr.1 == 1,
         B2Rd.Access == 1) %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "Accessible & low infection risk Scenario 1")

# 6. accumulation by abundance in low myrtle rust risk Scenario 2 and accessible areas
cumsum_list[[6]] <-
  cumsum_list[[6]] %>%
  filter(B2IRS2.lyr.1 == 1,
         B2Rd.Access == 1) %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "Accessible & low infection risk Scenario 2")

# combine
cumsum_df <- do.call(bind_rows, cumsum_list)

#Backup dataframe to file
write.csv(cumsum_df,"Model-outputs/Cumulative_area_scenarios.csv")

Abundance_area_summary <- cumsum_df %>%
  group_by(type) %>%
  summarise('Max area' = max(area_cumcum, na.rm = TRUE),
            'Max trees' = max(mean_cumsum, na.rm = TRUE))

write.csv(Abundance_area_summary[,1:3],"Model-outputs/Abundance_area_summary.csv")

#Re-order the factor levels for graphical display
cumsum_df$type = factor(cumsum_df$type,levels=c("All",
                                                "Accessible",
                                                "Low infection risk Scenario 1",
                                                "Low infection risk Scenario 2",
                                                "Accessible & low infection risk Scenario 1",
                                                "Accessible & low infection risk Scenario 2")) 

# plot - have updated x values of thresholds
cumuplot<- ggplot(data=subset(cumsum_df, type == "All" | 
                                type == "Low infection risk Scenario 1" |
                                type == "Low infection risk Scenario 2" | 
                                type == "Accessible")) +
            #geom_segment(aes(x = 2388, y = 0, xend = 2388, yend = 610), linetype = 4, color = "#000000",linewidth = 1) +
            #geom_segment(aes(x = 1244, y = 0, xend = 1244, yend = 220), linetype = 2, color = "#012840",linewidth = 1) +
            #geom_segment(aes(x = 1023, y = 0, xend = 1023, yend = 55), linetype = 3, color = "#025373",linewidth = 1) +
            #geom_segment(aes(x = 369, y = 0, xend = 369, yend = 11),linetype = 1,color="#03738C",linewidth = 1) +  
            #geom_segment(aes(x = 1023, y = 0, xend = 1023, yend = 55), linetype = 3, color = "#3FA8BF",linewidth = 1) +
            #geom_segment(aes(x = 369, y = 0, xend = 369, yend = 11),linetype = 1,color="#96D2D9",linewidth = 1) +   
            geom_line(aes(area_cumcum, mean_cumsum,
                  color = type), linewidth = 1.2) +
            #facet_wrap(~type,nrow=3) +
            labs(x = expression(paste("Cumulative land area (",km^2,")")),
                  y = "Cumulative abundance", color = "Management scenario") +
            coord_cartesian(expand = TRUE,
                  clip = "off",
                  ylim = c(0, NA),
                  xlim = c(0, 5222)) +
            scale_color_manual(values = c("#000000","#025373","#3FA8BF","#96D2D9")) +
            theme_classic() +
            theme(legend.position = c(0.75,0.3))

ggsave(cumuplot,
        filename=paste0("Model-outputs/Figure6.png"),
        device="png",
        height=10, width = 16, units = "cm",dpi="print")



#=== Distance decay in the Matern spatial covariance ===#
corplot <- plot(spde.posterior(spat_model, "shared_spatial", what = "matern.correlation"))
FigS7<- plot(corplot) +
          labs(x = "Distance (m)",
          y = "Residual correlation") +
          coord_cartesian(expand = FALSE,
                  ylim = c(0, 1)) +
          theme_classic()

ggsave(FigS7,
       filename=paste0("Model-outputs/FigureS7.png"),
       device="png",
       height=10, width = 16, units = "cm",dpi="print")