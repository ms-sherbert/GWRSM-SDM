# Modification of Hao Ran's inlabru script to explore the feasibility of 
# integrated SDM using point process model
# Written by: HR Lai and SM Herbert
# Written for R version 4.3.1 (SMH) and 4.2.2 (HRL)


#=== Preamble - use bits as needed for your computer ===#

rm(list=ls())
#Remember that working directory needs to be set to the local copy of the GWRSM-SDM repository

#Use if needed to install PointedSDMs package

#If running R v4.3.2 or more recent; presumably install from CRAN works
#install.packages('PointedSDMs')

#If running older version of R (e.g. 4.3.1)
#devtools::install_github("PhilipMostert/PointedSDMs")

#install.packages("INLA",repos=c(getOption("repos"),
#                INLA="https://inla.r-inla-download.org/R/stable"),
#                dep=TRUE)

#Load required packages
library(PointedSDMs) #version 1.3
require(inlabru)     #version 2.10.1
require(fmesher)
require(R6)
library(sf)
library(INLA)        #version 23.04.24 (HR Lai's computer) or version 23.09.09 (VUW PC)
library(terra)
library(tidyverse)
library(scoringRules)
library(tidyterra)


#=== Read in Data sources ===#

# Study domain
freshwater <-
  st_read("LCDB5-open-water/LCDB5-open-water-and-rivers.shp")
GWR <-
  st_read("GWRboundary/GWRboundary2193.shp") %>%
  st_difference(y = st_union(freshwater))

# Covariates
wet1 <- rast("GISinputs-repositories/NZEnvDS_v1.1/final_layers_nztm/topo_wetness.tif")
wet1 <- crop(wet1, GWR)
wet <- scale(wet1)

drain1 <- rast("GISinputs-repositories/NZEnvDS_v1.1/final_layers_nztm/soil_drainage.tif")
drain1 <- crop(drain1, GWR)
drain <- scale(drain1)

Tmin1 <- rast("GISinputs-repositories/NZEnvDS_v1.1/final_layers_nztm/temp_minColdMonth.tif")
Tmin1 <- crop(Tmin1, GWR)
Tmin <- scale(Tmin1)

precip1 <- rast("GISinputs-repositories/NZEnvDS_v1.1/final_layers_nztm/precip_warmQtr.tif")
precip1 <- crop(precip1, GWR)
precip <- scale(precip1)

humid1 <- rast("GISinputs-repositories/NZEnvDS_v1.1/final_layers_nztm/humidity_meanAnn.tif")
humid1 <- crop(humid1, GWR)
humid <- scale(humid1)

solar1 <- rast("GISinputs-repositories/NZEnvDS_v1.1/final_layers_nztm/solRad_winter.tif")
solar1 <- crop(solar1, GWR)
solar <- scale(solar1)

Trange1 <- rast("GISinputs-repositories/NZEnvDS_v1.1/final_layers_nztm/temp_annRange.tif")
Trange1 <- crop(Trange1, GWR)
Trange <- scale(Trange1)

road1 <- rast("GISinputs-repositories/NZEnvDS_v1.1/final_layers_nztm/distance_road.tif")
road1 <- crop(road1, GWR)
road <- scale(road1) #type in assigned name (e.g. 'road)' to check sf details

# SYZmai records

obs_all <-
  read_csv("SM_obs/SM-PO-PA.csv") %>%
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

#MRPM infection risk raster
IR <- rast("GISinputs-repositories/MR-IR-70/MRRisk70r.tif")
IR <- crop(IR, GWR)

#Distance to road raster
EUCroad <- rast("GISinputs-repositories/Euclid_distance_road/dis2road7Euc500r.tif")
EUCroad <- crop(EUCroad, GWR)

#=== Graph data sources for Fig2 of manuscript ===#

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
  cutoff = 500,   # need to revise later
  crs = st_crs(GWR)
)

mesh1$n
# plot(mesh1)
ggplot() +
  gg(data = GWR, fill = "darkgreen", colour = NA, alpha = 0.2) +
  gg(mesh1) +
  gg(data = subset(PAobs, NPres== 0), pch = 21, fill = "white") +
  gg(data = subset(PAobs, NPres == 1), pch = 21, fill = "black") +
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
covariates_minimal <-c(solar,drain,road) #minimal model containing only significant effects in global model


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

#Export model coefficients for fixed effects
write.csv(spat_model$summary.fixed,"Model-outputs/iPPM_fixed_effects.csv")

#=== Cross validation of model ===#
# By Leave-one-out cross validation (i.e. examining the effect of leaving out data types)
# This chunk currently doesn't work because the spat_model input is a copy model
# Needs fixing at intModel()
# Although, I'm not entirely sure that fixing this is a priority

#PA_out <- datasetOut(spat_model, dataset = 'PAobs', predictions = TRUE)


#==== Make predictions ====#
# The fixed factor covariate 'distance from road' is conditioned out of the predictive model (Warton et al. 2013)
# This is to remove the directional bias from our predictions
# i.e. make predictions over the whole range that would be true if every point in our space is equidistant from a road

projections <- predict(spat_model,
                       mesh=mesh1,
                       mask = GWR,
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

#=== Plot predictions for Figure 3c ===#
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
write.csv(scaled_df,"Model-outputs/iPPM_scaled_predictions.csv")
#Export projected values to a shapefile for analysis in ArcGIS or QGIS
st_write(scaled_df,"Model-outputs/iPPM_predictions.shp",append=FALSE)

#Calculate how many scaled points were predicted to be greater than the average 
length(which(scaled_df$scaled_mean > 0))

#Aggregate predictions to 1km grid and export to shapefile
#New_grid <-
#  st_read("GISinputs-repositories/GWR_1km_MRPM_grid/GWR_1km_MRPM_grid.shp")
#pred_agg <- aggregate(projections$predictions,New_grid,FUN = sum)
#st_write(pred_agg,"Model-outputs/Agg_iPPM_predictions.shp")

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

write.csv(Lambda_total$predictions,"Model-outputs/iPPM_totabundance_predictions.csv")

# Generate posterior abundance distribution
# The 95% CrI of Lambda_total is used to define a broader range of plugin values
# to constrain posterior values of N
# We could also use a value just lower than the number of 268 observations to constrain the lower tail
# Because we know we won't have fewer than this number of trees
Trees <- predict(
  spat_model,
  fm_int(mesh1, GWR),
  formula = ~ data.frame(
    N = seq(100, 1000, 10),
    dpois(seq(100, 1000, 10),
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


#=== Prediction score ===#
# spat_model$bru_info$model$formula
# spat_model$componentsJoint
# spat_model$bru_info$lhoods$POobs_geometry$inla.family
# spat_model$bru_info$lhoods$PAobs_NPres$inla.family

# Integrated predictions
# First we partition the landscape into grids of some resolution
source("Code/utils.R")
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

B2IR <- terra::extract(IR, B2Pts)
#Change values to binary multipliers (0 = unacceptably high risk, 1 = acceptable risk)
B2IR["MRRisk70r"][B2IR["MRRisk70r"] == 100] <- 1
B2IR["MRRisk70r"][B2IR["MRRisk70r"] == 128] <- NA

B2Rd <- terra::extract(EUCroad, B2Pts)
Access <- as.numeric(B2Rd$Area)
B2Rd <- cbind(B2Rd,Access)
#Change values to binary multipliers (0 = inaccessible, 1 = accessible)
B2Rd["Access"][B2Rd["Access"] == 1] <- 0
B2Rd["Access"][B2Rd["Access"] == 2] <- 1

B2 <-cbind(B2,B2IR$MRRisk70r,B2Rd$Access)

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

# construct cumulative curve in three ways
cumsum_list <- list(abun_total, abun_total, abun_total)

# 1. accumulation by abundance only
cumsum_list[[1]] <-
  cumsum_list[[1]] %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "All suitably waterlogged soils")

# 2. accumulation by both abundance in low myrtle rust risk areas only
cumsum_list[[2]] <-
  cumsum_list[[2]] %>%
  filter(B2IR.MRRisk70r == 1) %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "Infection risk < 0.7")

# 3. accumulation by both abundance in accessible areas only
cumsum_list[[3]] <-
  cumsum_list[[3]] %>%
  filter(B2Rd.Access == 1) %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "Accessible")

# 4. accumulation by both abundance in low myrtle rust risk and accessible areas
cumsum_list[[4]] <-
  cumsum_list[[4]] %>%
  filter(B2IR.MRRisk70r == 1,
         B2Rd.Access == 1) %>%
  # sort by mean abundance per block in decreasing order
  arrange(desc(mean)) %>%
  # calculate cumulative sum
  mutate(mean_cumsum = cumsum(mean),
         area = as.numeric(st_area(.))/1e6,
         area_cumcum = cumsum(area),
         type = "Accessible & infection risk < 0.7")

# combine
cumsum_df <- do.call(bind_rows, cumsum_list)

#Backup dataframe to file
write.csv(cumsum_df,"Model-outputs/Cumulative_area_scenarios.csv")

#Re-order the factor levels for graphical display
cumsum_df$type = factor(cumsum_df$type,levels=c("All suitably waterlogged soils",
                                                "Infection risk < 0.7",
                                                "Accessible",
                                                "Accessible & infection risk < 0.7")) 

# plot
ggplot(cumsum_df) +
  geom_segment(aes(x = 2479, y = 0, xend = 2479, yend = 470), linetype = 4, color = "#000508",size = 1) +
  geom_segment(aes(x = 857, y = 0, xend = 857, yend = 230), linetype = 2, color = "#006d9b",size = 1) +
  geom_segment(aes(x = 1090, y = 0, xend = 1090, yend = 55), linetype = 3, color = "#00b3fc",size = 1) +
  geom_segment(aes(x = 390, y = 0, xend = 390, yend = 15),linetype = 1,color="#5fd0ff",size = 1) +  
  geom_line(aes(area_cumcum, mean_cumsum,
                color = type), size = 1.2) +
  labs(x = "Cumulative land area (sq. km)",
       y = "Cumulative abundance",
       color = "Management scenario") +
  coord_cartesian(expand = TRUE,
                  clip = "off",
                  ylim = c(0, NA),
                  xlim = c(0, 8500)) +
  scale_color_manual(values = c("#000508","#006d9b","#00b3fc","#5fd0ff")) +
  theme_bw() +
  theme(legend.position = c(0.81, 0.3))


#=== Distance decay in the Matern spatial covariance ===#
corplot <- plot(spde.posterior(spat_model, "shared_spatial", what = "matern.correlation"))
plot(corplot) +
  labs(x = "Distance (m)",
       y = "Residual correlation") +
  coord_cartesian(expand = FALSE,
                  ylim = c(0, 1)) +
  theme_classic()


#=== Create predictions of abundance for every 100x100 m ===#

# Note that this chunk won't run on a regular PC or laptop
# In theory it should run on a HPC, we are looking into this
# Partition the landscape into grids of 100 x 100 m resolution
# note that resolution defines the area of the grid squares in square metres
B3 <- partition(samplers = GWR, resolution = 100)
#plot(B2)
#st_write(B2,"GISinputs-repositories/B2_sampler.shp",append=FALSE)

# add the total observed number of occurrence points
B3$NPres <- lengths(st_intersects(B3, filter(obs_all, NPres != 0)))

Lambda <- predict(
  spat_model,
  fm_int(domain = mesh1, sampler = B3),  # this is to integrate the prediction over 1 x 1 km grid area
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


B3SpatVec <- vect(as(B3, "Spatial"))
B3Pts <- centroids(B3SpatVec, inside=FALSE) #currently a workaround to account for imperfect alignment of the grids

B3IR <- terra::extract(IR, B3Pts)
#Change values to binary multipliers (0 = unacceptably high risk, 1 = acceptable risk)
B3IR["MRRisk70r"][B3IR["MRRisk70r"] == 100] <- 1
B3IR["MRRisk70r"][B3IR["MRRisk70r"] == 128] <- NA

B3Rd <- terra::extract(EUCroad, B3Pts)
Access <- as.numeric(B3Rd$Area)
B3Rd <- cbind(B3Rd,Access)
#Change values to binary multipliers (0 = inaccessible, 1 = accessible)
B3Rd["Access"][B3Rd["Access"] == 1] <- 0
B3Rd["Access"][B3Rd["Access"] == 2] <- 1

B3 <-cbind(B3,B3IR$MRRisk70r,B3Rd$Access)

# join B3 to predicted abundances
abun_total <-
  B3 %>%
  rownames_to_column(".block") %>%
  left_join(
    Lambda$predictions %>%
      rownames_to_column(".block")
  )

# Scale and write the abundance per 100 x 100 km cell to shapefile
scaled_mean_100m <- scale(abun_total$mean)
abun_total<-cbind(abun_total,scaled_mean_1km)
st_write(abun_total,"Model-outputs/Abundance100m.shp") 
plot(abun_total[20]) #