# GWRSM-SDM

Repository for the manuscript Herbert SM*, Tomscha SA*, Lai HR, Benavidez R, Balkwill CG, Ruston P, Jackson B, Deslippe JR (Accepted April 2025) Identifying potentially suitable and accessible refugia to mitigate impacts of an emerging disease on a rare tree. Conservation Biology.
Authors labelled with an asterisk (*) contributed equally. 

Please note that this manuscript is currently in press. 

[![DOI](https://zenodo.org/badge/772350956.svg)](https://zenodo.org/doi/10.5281/zenodo.10836317)

## Repository contents

`Code` folder contains all R code associated with this project.

`GWRboundary` contains the shapefile of the Greater Wellington Regional boundary. 

`LCDB5-open-water` contains a shapefile of open water (ponds, lakes and major rivers) in the Greater Wellington Region. This is derived from the Landcover Database version 5 (LCDB5), merged from the 'Lakes and Ponds' and 'Rivers' 2018 cover classes. Projection: NZTM2000 (EPSG 2913). 

`Model-outputs` contains outputs from the most recent runs of the codes for the point process models. 

`SM_obs_public` contains presence-only and presence-absence observation data for swamp maire. The file `SM-PO-PA.csv` contains all the swamp maire data and `SMobs269.csv` contains presence-only swamp maire data. The file `GWRmyrtles.csv` contains presence-only swamp maire data plus records of an additional 19 myrtle species for comparison in bias assessment. 

## Useage

### Dependencies

All code provided in this repository has been tested for R versions 4.3.1 and 4.2.2. The versions of key R packages used were `PointedSDMs` version 1.3, `inlabru` version 2.10.1, `INLA` versions 23.04.24 and 23.09.09, `terra` versions 1.7-46 (most recent tests) and 1.7-55, `sf` version 1.0-14, and `sampbias` version 1.0.5. We cannot guarantee that these codes will work with other versions of R or these packages.  

### Code files

To produce the analyses reported in the manuscript, we ran the codes in the following order:

1. `AssessSpatBias.R` and `QuantSampBias.R` were used to assess bias in the presence-only swamp maire data.
2. `PPM_explore.R` and `PointedSDMS_Waitrial.R` were used for presence-only point process modeling and integrated point process modeling, respectively. Both codes have a dependency on `utils.R`.
3. `Bootstrap_LUCI_filter.R` was used to determine the habitat suitability of different LUCI flood risk classifications for swamp maire, and to filter the integrated point process model predictions.
4. `Graph-model-effects.R` uses the `Model-outputs/Output-summaries/Model_coefficients.csv` file to plot point process model covariate effects for Figure 2a of the manuscript (see: 10.6084/m9.figshare.25404670). 
5. `Refugia-plot.R` uses the `Model-outputs/Output-summaries/Restorable_areas.csv` output from the integrated point process models to create Figure 5 of the manuscript (see: 10.6084/m9.figshare.25405432). 

### Datasets called in by codes

Many of the data files used to run the codes have been provided in full in this repository. However, there are some important exceptions as follows:

We have omitted geospatial coordinates from the data rows that contain sensitive location data in the folder `SM_obs_public`. Omitted rows can be recognized by having a latitude of 170 and longitude of -40. Please contact the corresponding author julie.deslippe[at]vuw.ac.nz to request access to the full data set so that we can obtain the necessary permissions. 

Large geospatial input files are stored separately. These are .tif raster files and can be downloaded from FigShare:

- floodrisk2023.tif: 10.6084/m9.figshare.25413709
- Scenario1.tif and Scenario2.tif: 10.6084/m9.figshare.25413733
- dis2road7Euc500r.tif: 10.6084/m9.figshare.25413748

The raster files of bioclimatic covariates used in the point process model codes (`PPM_explore.R` and `PointedSDMs_Waitrial.R`) are from the New Zealand Environmental Data Stack (NZEnvDS) version 1.1 (NZTM projection), and can be downloaded here: https://doi.org/10.7931/m6rm-vz40. Should you use these data, please cite the following paper and any others of relevance:

> McCarthy, J. K., Leathwick, J. R., Roudier, P., Barringer, J. R. F., Etherington, T. R., Morgan, F. J., Odgers, N. P., Price, R. H., Wiser, S. K., Richardson, S. J. (2021) New Zealand Environmental Data Stack (NZEnvDS): A standardised collection of environmental spatial layers for biodiversity modelling and site characterisation. New Zealand Journal of Ecology, 45(2): 3440 https://newzealandecology.org/nzje/3440


