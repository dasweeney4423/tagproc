polyprox <- function(data, polypath) {
  # Format pos/bhvr data
  LAT <- data[, which(names(data) %in% c("Latitude", "StartLat"))]
  LON <- data[, which(names(data) %in% c("Longitude", "StartLon"))]

  polylf <- list.files(polypath)
  if (length(polylf) == 0) {
    polyname <- polypath
    # Load polygon
    polj <- suppressWarnings(rgdal::readOGR(dsn = polyname))
    polyname <- strsplit(polyname,'/')[[1]][length(strsplit(polyname,'/')[[1]])]

    # Loop through each GPS record
    polydist <- rep(NA, times = nrow(data))
    for (i in 1:nrow(data)) {
      # Collect lat and long
      x <- LON[i]
      y <- LAT[i]

      if (!any(is.na(c(x, y)))) {
        # Format spatial datasets
        coord <- data.frame(longitude = x, latitude = y)
        sp::coordinates(coord) <- c("longitude", "latitude")
        sp::proj4string(coord) <- sp::proj4string(polj)

        # Test whether point is in polygon (if FALSE, it is not)
        inpoly <- !is.na(sp::over(coord, as(polj, "SpatialPolygons")))
        if (inpoly) {
          polydist[i] <- 0 # If pt is in polygon, set distance to 0
        } else {
          dij <- geosphere::dist2Line(p = coord, line = polj)  # If pt is not in polygon, determine distance to polygon
          distij <- round(dij[1] / 1000, digits = 3)  # Format distance
          polydist[i] <- distij
        }
      } # end of if x | y is na
    } # end of GPS loop

    # Add distance vector to GPS dataframe
    data$poly <- polydist
    names(data)[which(names(data) == "poly")] <- polyname
  } else {
    # Loop through each Polygon
    for (j in 1:length(polylf)) {
      # Load polygon j
      polyname <- paste0(polypath, "/", polylf[j])
      polj <- suppressWarnings(rgdal::readOGR(dsn = polyname))
      polyname <- strsplit(polyname,'/')[[1]][length(strsplit(polyname,'/')[[1]])]

      # Loop through each GPS record
      polydist <- rep(NA, times = nrow(data))
      for (i in 1:nrow(data)) {
        # Collect lat and long
        x <- LON[i]
        y <- LAT[i]

        if (!any(is.na(c(x, y)))) {
          # Format spatial datasets
          coord <- data.frame(longitude = x, latitude = y)
          sp::coordinates(coord) <- c("longitude", "latitude")
          sp::proj4string(coord) <- sp::proj4string(polj)

          # Test whether point is in polygon (if FALSE, it is not)
          inpoly <- !is.na(sp::over(coord, as(polj, "SpatialPolygons")))
          if (inpoly) {
            polydist[i] <- 0 # If pt is in polygon, set distance to 0
          } else {
            dij <- geosphere::dist2Line(p = coord, line = polj)  # If pt is not in polygon, determine distance to polygon
            distij <- round(dij[1] / 1000, digits = 3)  # Format distance
            polydist[i] <- distij
          }
        } # end of if x | y is na
      } # end of GPS loop

      # Add distance vector to GPS dataframe
      data$poly <- polydist
      names(data)[which(names(data) == "poly")] <- polyname
    } # end of polygon loop
  }

  return(data)
}
