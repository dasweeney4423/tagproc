

freitas.filter <- function(argosfile, deployfile = NULL,
                           vmax = 2.778, ang = c(15, 25), distlim = c(3000, 3000)) {
  #load argos file and deployment data if desired
  locations <- load.argos(argosfile, deployfile)
  locations <- locations[!is.na(locations$Latitude) & !is.na(locations$Longitude),]

  if (nrow(locations) > 0) {
    D <- as.character(locations$Date)

    filter <- argosfilter::sdafilter(lat = locations$Latitude,
                                        lon = locations$Longitude,
                                        dtime = D,
                                        lc = locations$LocationQuality,
                                        vmax = vmax,
                                        ang = ang,
                                        distlim = distlim)

    locations <- locations[filter != "removed",]

    # Add column with info on source file
    locations$argos.source <- rep(source.file, times = nrow(locations))
  }

  return(locations)
}
