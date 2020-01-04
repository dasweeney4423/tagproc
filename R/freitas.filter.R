freitas.filter <- function(argos, vmax, ang = c(15, 25), distlim = c(3000,3000)) {
  argos <- argos[!is.na(argos$Latitude) & !is.na(argos$Longitude),]

  if (nrow(argos) > 0) {
    D <- as.character(argos$Date)
    filter <- argosfilter::sdafilter(lat = argos$Latitude,
                                     lon = argos$Longitude,
                                     dtime = D,
                                     lc = argos$LocationQuality,
                                     vmax = vmax,
                                     ang = ang,
                                     distlim = distlim)
    argos <- argos[filter != "removed",]
  }

  return(argos)
}
