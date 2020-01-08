#' Perform a simple GPS filter
#'
#' This function filters GPS location data using a specified threshold for the number of satellites and residual values
#' @param gps A dataframe of gps locations to be sent through the filter
#' @param min.sat Minimum number of satellites for that location for it to be kept. Default is 5
#' @param max.residual Maximum residual value. For a location to be kept, its residual must be less than or equal to this input. Default is 35 (Dujon et al. 2014)
#' @return A dataframe of filtered locations that no longer includes outliers.
#' @examples #examples not yet provided, sorry :(

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
