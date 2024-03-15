partition <- function(samplers, resolution = NULL, nrows = NULL, ncols = NULL) {
  # Create a grid for the given boundary
  if (is.null(resolution)) {
    grid <- rast(terra::ext(samplers),
                 crs = st_crs(samplers),
                 nrows = nrows, ncols = ncols
    )
  }

  if (is.null(c(nrows, ncols))) {
    grid <- rast(terra::ext(samplers),
                 crs = st_crs(samplers),
                 resolution = resolution
    )
  }

  gridPolygon <- terra::as.polygons(grid)

  # Extract the boundary with subpolygons only
  out <-
    sf::st_as_sf(terra::intersect(gridPolygon, terra::vect(samplers))) %>%
    sf::st_set_crs(st_crs(samplers))
}

