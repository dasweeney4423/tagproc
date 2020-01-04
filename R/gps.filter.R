gps.filter <- function(gps, min.sat = 5, max.residual = 35) {
  gps <- gps[!is.na(gps$Latitude) & !is.na(gps$Longitude),]

  if (nrow(gps) > 0) {
    bads <- which(!is.na(gps$Bad.Sats))
    gps$Good.Satellites <- gps$Satellites
    gps$Good.Satellites[bads] <- gps$Satellites[bads] - gps$Bad.Sats[bads]
    gps <- gps[gps$Good.Satellites >= min.sat,]
    res.test <- which(gps$Residual <= max.residual)
    if (length(res.test) > 0) {gps <- gps[res.test,]}
  }

  return(gps)
}
