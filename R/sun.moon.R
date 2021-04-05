#' Sync environmental data with tag data
#'
#' Sync environmental data with tag location and time data from the WC portal
#' @param data A dataframe of tag locations and times to which environmental data will be attached
#' @return A dataframe with all the same columns as the input data but extra columns containing enironemtnal data are attached.
#' @examples #examples not yet provided, sorry :(

sun.moon <- function(data) {
  data <- data.frame(data)
  TIME <- time.turner(data[, which(names(data) %in% c("Date", "StartTime", 'Day'))])$strp
  LAT <- data[, which(names(data) %in% c("Latitude", "StartLat", "Lat"))]
  LON <- data[, which(names(data) %in% c("Longitude", "StartLon", "Long"))]

  # Setup variables
  sunrise <- sunset <- solarNoon <- civilDawn <- endDawn <- civilDusk <- startDusk <- rep(NA, times = nrow(data))
  sunAzimuth <- sunAltitude <- moonAzimuth <- moonAltitude <- moonIlluminatedFraction <- moonPhase <- rep(NA, times = nrow(data))

  # Loop thru each entry, calculate variables
  for (i in 1:nrow(data)) {
    posmatrix <- matrix(c(LON[i], LAT[i]), nrow = 1) # Create GPS matrix
    timei <- as.POSIXct(format(TIME[i], tz = "UTC"), tz = "UTC")

    if ((!is.na(posmatrix[1])) & (!is.na(posmatrix[2]))) {
      sunrise[i] <- as.character(as.POSIXct(maptools::sunriset(posmatrix, timei, direction = "sunrise", POSIXct.out = T)$time))
      sunset[i]  <- as.character(as.POSIXct(maptools::sunriset(posmatrix, timei, direction = "sunset", POSIXct.out = T)$time))
      solarNoon[i] <- as.character(as.POSIXct(maptools::solarnoon(posmatrix, timei, POSIXct.out = T)$time))
      civilDawn[i] <- as.character(as.POSIXct(maptools::crepuscule(posmatrix, timei, solarDep = 6, direction = "dawn", POSIXct.out = T)$time))
      endDawn[i] <- as.character(as.POSIXct(maptools::crepuscule(posmatrix, timei, solarDep = -6, direction = "dawn", POSIXct.out = T)$time))
      civilDusk[i] <- as.character(as.POSIXct(maptools::crepuscule(posmatrix, timei, solarDep = 6, direction = "dusk", POSIXct.out = T)$time))
      startDusk[i] <- as.character(as.POSIXct(maptools::crepuscule(posmatrix, timei, solarDep = -6, direction = "dusk", POSIXct.out = T)$time))

      sunangle <- oce::sunAngle(timei, posmatrix[1, 1], posmatrix[1, 2])  #package oce and maptools packages calculate the sun and moon algortihmically. No lookup files needed.
      sunAzimuth[i] <- sunangle$azimuth
      sunAltitude[i] <- sunangle$altitude
      moonangle <- oce::moonAngle(timei, posmatrix[1, 1], posmatrix[1, 2])
      moonAzimuth[i] <- moonangle$azimuth
      moonAltitude[i] <- moonangle$altitude
      moonIlluminatedFraction[i] <- moonangle$illuminatedFraction
      moonPhase[i] <- moonangle$phase %% 1
    }
  }

  # Add calculated variables to results dataframe
  data$sunrise <- time.turner(sunrise)$raw
  data$sunset <- time.turner(sunset)$raw
  data$solarNoon <- time.turner(solarNoon)$raw
  data$civilDawn <- time.turner(civilDawn)$raw
  data$endDawn <- time.turner(endDawn)$raw
  data$civilDusk <- time.turner(civilDusk)$raw
  data$startDusk <- time.turner(startDusk)$raw
  data$sunAzimuth <- sunAzimuth
  data$sunAltitude <- sunAltitude
  data$moonAzimuth <- moonAzimuth
  data$moonAltitude <- moonAltitude
  data$moonIlluminatedFraction <- moonIlluminatedFraction
  data$moonPhase <- moonPhase

  return(data)
}
