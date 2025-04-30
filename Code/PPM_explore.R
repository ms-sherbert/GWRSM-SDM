# Script to explore the feasibility of point process model
# Written by: HR Lai & S Herbert
# Written for R version 4.2.2 (HR Lai), but works fine with 4.3.1 (S Herbert)

#=== Preamble - use bits as needed for your computer ===#

rm(list=ls())

local.dir <- "D:/" #Change so that file path matches where these files are saved on your computer
local.files <- "D:/Repositories/Offline-files-GWRSM-SDM/" #file path where any gitignored files are stored locally

#Remember that working directory needs to be set to the local copy of the GWRSM-SDM repository, e.g.
setwd(paste0(local.dir,"Repositories/GWRSM-SDM"))

run.date <- as.character(Sys.Date())

#=== Load required packages ===#

library(sf)
library(terra)
library(tidyverse)
library(inlabru) #version 2.10.1
library(INLA) #version 23.04.24 (HR Lai's computer) or version 23.09.09 (VUW PC)
library(scoringRules)
library(tidyterra)
source("Code/utils.R")

#=== Read in Data sources ===#

# Study domain
freshwater <-
  st_read("LCDB5-clipping-layers/LCDB5-open-water-and-rivers.shp")
GWR <-
  st_read("GWRboundary/GWRboundary2193.shp") %>%
  st_difference(y = st_union(freshwater)) 

# Note that repository file paths are currently set to work from VUW PC
# File paths for Hao Ran's computer follow this format: "covs_SM_Feb2023/distance_road.tif"

wet <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/topo_wetness.tif"))
wet <- crop(wet, GWR)
wet <- scale(wet)

drain <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/soil_drainage.tif"))
drain <- crop(drain, GWR)
drain <- scale(drain)

Tmin <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/temp_minColdMonth.tif"))
Tmin <- crop(Tmin, GWR)
Tmin <- scale(Tmin)

precip <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/precip_warmQtr.tif"))
precip <- crop(precip, GWR)
precip <- scale(precip)

humid <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/humidity_meanAnn.tif"))
humid <- crop(humid, GWR)
humid <- scale(humid)

solar <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/solRad_winter.tif"))
solar <- crop(solar, GWR)
solar <- scale(solar)

Trange <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/temp_annRange.tif"))
Trange <- crop(Trange, GWR)
Trange <- scale(Trange)

road <- rast(paste0(local.files,"NZEnvDS_v1.1/final_layers_nztm/distance_road.tif"))
road <- crop(road, GWR)
road <- scale(road) #type in assigned name (e.g. 'road)' to check sf details

# res(Tmin); res(precip); res(humid); res(solar); res(frost); res(Trange); res(road) #to check resolution of rasters


#=== SYZmai records ===#

# Presence-only records (pooled from NVS, iNaturalist, authors' observations, and herbaria)
obs <-
  read_csv(paste0(local.files,"SMobs269.csv")) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = "+proj=longlat +ellips=WGS84") %>%
  st_transform(crs = st_crs(GWR)) %>%
  st_intersection(y = GWR)


#=== Model preparation ===#

# Mesh
# following https://rpubs.com/jafet089/886687
# also https://haakonbakkagit.github.io/btopic104.html
coords <- st_coordinates(obs)
domain <- as(GWR, "Spatial") %>% inla.sp2segment()
# domain <- as_Spatial(GWR)
# max.range <- max(apply(coords, 2, function(x) diff(range(x))))
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

mesh1$n #Number of vertices in mesh
# plot(mesh1)
ggplot() +
  # gg(data = flood) +
  # gg(data = Tmin) +
  gg(data = GWR, fill = "darkgreen", colour = NA, alpha = 0.2) +
  gg(mesh1) +
  gg(data = obs, pch = 21, fill = "purple") +
  ggtitle(label = "A") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

saveRDS(mesh1,paste0("Model-outputs/poPPM_mesh_GWR_",run.date,".rds"))
#mesh1 <- read_rds("Model-outputs/poPPM_mesh_GWR_20240815.rds")

#================= Modelling ==================#

# spatial covariance component
# priors are the same as the integrated PPM
matern <- inla.spde2.pcmatern(
  mesh1,
  prior.sigma = c(1, 0.5),
  prior.range = c(1000, 0.25)
)

# combine all model components together
cmp <- geometry ~
  fmatern(geometry, model = matern) +
  froad(road, model = "linear") +
  fTmin(Tmin, model = "linear") +
  fprecip(precip, model = "linear") +
  fhumid(humid, model = "linear") +
  fsolar(solar, model = "linear") +
  fTrange(Trange, model = "linear") +
  fwet(wet, model = "linear") +
  fdrain(drain, model = "linear") +
  Intercept(1)


fit <- lgcp(
  cmp,
  data = obs,
  domain = list(geometry = mesh1),
  options = list(control.inla = list(int.strategy='eb'))
)

saveRDS(fit,paste0("Model-outputs/poPPM_GWR_",run.date,".rds"))

# code to read from saved model if needed
#fit <- read_rds("Model-outputs/poPPM_GWR_2024-08-19.rds")

summary(fit)

#Export model coefficients for fixed effects
write.csv(fit$summary.fixed,paste0("Model-outputs/poPPM_fixed_effects_",run.date,".csv"))

#=== Calculate prediction score ===#
# Integrated predictions
# First we partition the landscape into grids of some resolution (in km * km)

B1 <- partition(samplers = GWR, resolution = 50*50)
# plot(B1)

# add the total observed number of occurrence points
B1$NPres <- lengths(st_intersects(B1, obs))

# predict the total number of occurrences (counts) per grid
# As per our prediction of total abundance in the GWR,
# We leave the following line out fo this chunk:
# PAobs_intercept +
# This is so that the model only predicts the number of points where the species is present

Lambda <- predict(
  fit,
  fm_int(domain = mesh1, sampler = B1),  # this is to integrate the prediction over grid area
  formula = ~ tapply(weight * exp(fmatern +
                                    fTmin +
                                    fprecip +
                                    fhumid +
                                    fsolar +
                                    fwet +
                                    fdrain +
                                    fTrange +
                                    froad + #For CRPS we want road back in
                                    Intercept),
                     .block,
                     sum)
)

abun_total <-
  Lambda %>%
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

#=== Examine Parameters ====#
# int.plot <- plot(fit, "Intercept")
flist <- vector("list", NROW(fit$summary.random$fTmin))
for (i in seq_along(flist)) flist[[i]] <- plot(fit, "fTmin", index = i)
multiplot(plotlist = flist, cols = 3)

spde.range <- spde.posterior(fit, "fmatern", what = "range")
spde.logvar <- spde.posterior(fit, "fmatern", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)
multiplot(range.plot, var.plot)

corplot <- plot(spde.posterior(fit, "fmatern", what = "matern.correlation"))
covplot <- plot(spde.posterior(fit, "fmatern", what = "matern.covariance"))
multiplot(covplot, corplot)

#==== Make predictions ====#

pred <- predict(
  fit,
  #fm_pixels(mesh1, mask = GWR),
  fm_pixels(mesh=mesh1, mask = GWR, dims = c(1445,1021)), 
  ~ data.frame(
    lambda = exp(fmatern + fTmin + fprecip + fhumid + fsolar + fwet +fdrain + fTrange + Intercept),
    loglambda = fmatern + fTmin + fprecip + fhumid + fsolar + fwet + fdrain + fTrange + Intercept,
    log_spatres = fmatern
  ),
  num.threads = NULL
)

write_rds(pred, paste0("Model-outputs/poPPM_predictions_",run.date,".rds"))
#pred <- read_rds("Model-outputs/poPPM_predictions_2024-08-19.rds")

#=== Plot predictions for Figure 2 ===#

pl1 <- ggplot() +
  gg(pred$lambda, geom = "tile") +
  geom_sf(data = GWR, alpha = 0.1) +
  #geom_sf(data = obs) +
  scale_fill_viridis_c(name = "Predicted N records") +
  ggtitle("LGCP fit to Points", subtitle = "(Response Scale)") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

#I edited pl2 a bit to create FigS5a
pl2 <- ggplot() +
  gg(pred$loglambda, geom = "tile",alpha = 1.0) +
  geom_sf(data = GWR, alpha = 0.0) +
  #geom_sf(data = obs) +
  labs(xlab="",ylab="")+
  scale_fill_viridis_c(name = "Predicted mean intensity") +
  ggtitle("A") +
  xlab("Longitude") +
  ylab("Latitude") +
  #ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)")+
  theme_bw()

pl3 <- ggplot() +
  gg(pred$log_spatres, geom = "tile") +
  geom_sf(data = GWR, alpha = 0.1) +
  geom_sf(data = obs) +
  scale_fill_viridis_c() +
  ggtitle("LGCP fit to Points", subtitle = "(Matern residuals)")+
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

# multiplot(pl1, pl2, cols = 2)
# multiplot(pl2, pl3, cols = 2)
multiplot(pl1, pl3, pl2, cols = 2)


#scaled means plot for Fig 3b
scaled_mean_PO <- scale(pred$loglambda$mean)
scaled_df_PO <-cbind(scaled_mean_PO,pred$loglambda)

ScaledPreds_PO <- ggplot() +
  geom_sf(data=scaled_df_PO, aes(geometry = geometry, color = scaled_mean_PO)) +
  geom_sf(data = GWR, alpha = 0.0) +
  #geom_sf(data = subset(PAobs, NPres == 0), color = "black") +
  #geom_sf(data = POobs, color = "grey70",pch=15) +
  #geom_sf(data = subset(PAobs, NPres == 1), color = "grey70") +
  ggtitle("B") +
  scale_color_viridis_c(name = "Predicted mean intensity (scaled)") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

#plot of standard deviation

StDevs_PO <- ggplot() +
  geom_sf(data=scaled_df_PO, aes(geometry = geometry, color = sd)) +
  geom_sf(data = GWR, alpha = 0.0) +
  #geom_sf(data = subset(PAobs, NPres == 0), color = "black") +
  #geom_sf(data = POobs, color = "grey70",pch=15) +
  #geom_sf(data = subset(PAobs, NPres == 1), color = "grey70") +
  ggtitle("Standard deviation of predicted mean intensity - presence-only PPM") +
  scale_color_viridis_c(name = "Standard deviation") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

#Export projected values to a csv file
write.csv(scaled_df_PO,paste0("Model-outputs/PO-PPM_scaled_predictions_100_",run.date,".csv"))
#Export projected values to a shapefile for analysis in ArcGIS or QGIS
st_write(scaled_df_PO,paste0("Model-outputs/PO-PPM_predictions_100_",run.date,".shp"))

ggsave(pl2,
       filename=paste0("Model-outputs/PO-LinearPreds100m.png"),
       device="png",
       height=10, width = 16, units = "cm",dpi="print")

ggsave(ScaledPreds_PO,
       filename=paste0("Model-outputs/PO-ScaledPreds100m.png"),
       device="png",
       height=10, width = 16, units = "cm",dpi="print")

ggsave(StDevs_PO,
       filename=paste0("Model-outputs/PO-SdPreds100m.png"),
       device="png",
       height=10, width = 16, units = "cm",dpi="print")

#Calculate how many scaled points were predicted to be greater than the average 
length(which(scaled_df_PO$scaled_mean_PO > 0))

#=== Calculate total swamp maire abundance in the GWR ===#
# We didn't end up using this code for the paper, but have left it in for reference 
# See also https://inlabru-org.github.io/inlabru/articles/2d_lgcp_sf.html#estimating-abundance
# We leave the following line out of this code chunk:
# PAobs_intercept +
# This is so that the model only predicts the number of points where the species is present

Lambda_total <- predict(
  fit,
  fm_int(mesh1, GWR),
  formula = ~ sum(weight * exp(
    fmatern +
      fTmin +
      froad +
      fprecip +
      fhumid +
      fsolar +
      fwet +
      fdrain +
      fTrange +
      Intercept))
)

Lambda_total

write.csv(Lambda_total$predictions,"Model-outputs/PO-PPM_totabundance_predictions.csv")

# Generate posterior abundance distribution
# The 95% CrI of Lambda_total is used to define a broader range of plugin values
# to constrain posterior values of N
Trees <- predict(
  fit,
  fm_int(mesh1, GWR),
  ~ data.frame(
    N = 10:5200,
    dpois(10:5200,
          lambda = sum(weight * exp(
            fmatern +
              fTmin +
              fprecip +
              fhumid +
              fsolar +
              fwet +
              fdrain +
              fTrange +
              Intercept))
    )
  )
)

#Get quantiles of the posterior abundance distribution
inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x = Trees$N, y = Trees$mean))

#Get mean of the posterior abundance distribution
inla.emarginal(identity, marginal = list(x = Trees$N, y = Trees$mean))

#Plot posterior and plugin abundance distributions for comparison
Trees$plugin_estimate <- dpois(Trees$N, lambda = Lambda_total$mean)

PostN <- ggplot(data = Trees) +
  geom_line(aes(x = N, y = mean, colour = "Posterior")) +
  geom_line(aes(x = N, y = plugin_estimate, colour = "Plugin")) +
  theme_bw()