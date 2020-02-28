#' Perform a simple GPS filter
#'
#' This function filters GPS location data using a specified threshold for the number of satellites and residual values
#' @param gps A dataframe of gps locations to be sent through the filter
#' @param max.te Maximum time error (absolute value) desired kept. Default is 3
#' @param max.residual Maximum residual value. For a location to be kept, its residual must be less than or equal to this input. Default is 35 (Dujon et al. 2014)
#' @param min.sat Minimum number of satellites for that location for it to be kept. Default is no satellite filtering
#' @return A list of three elements in similar format to the data provided in the input for gps:
#' \itemize{
#' \item{\strong{all: }} A dataframe containing all input gps locations with an additional column specifying whether the filter marked the given location as an outlier or not.
#' \item{\strong{retained: }} A dataframe containing all retained input gps locations following the execution of the filter algorithm.
#' \item{\strong{outliers: }} A dataframe containing all filtered-out input gps locations following the execution of the filter algorithm.
#' }
#' @examples #examples not yet provided, sorry :(

gps.filter <- function(gps, max.te = 3, max.residual = 35, min.sat = NULL) {
  gps <- gps[!is.na(gps$Latitude) & !is.na(gps$Longitude),]

  if (nrow(gps) > 0) {
    bads <- which(!is.na(gps$Bad.Sats))
    gps$Good.Satellites <- gps$Satellites
    gps$Good.Satellites[bads] <- gps$Satellites[bads] - gps$Bad.Sats[bads]
    test <- gps
    if (!is.null(max.te)) {
      test <- test[which(abs(test$Time.Error) <= max.te),]
    }
    if (!is.null(max.residual)) {
      test <- test[which(test$Residual <= max.residual),]
    }
    if (!is.null(min.sat)) {
      test <- test[which(test$Good.Satellites >= min.sat),]
    }
  }

  for (i in 1:nrow(gps)) {
    if (gps$Date[i] %in% test$Date) {
      gps$outlier[i] <- FALSE
    } else {
      gps$outlier[i] <- TRUE
    }
  }
  retained <- gps[which(gps$outlier == FALSE),]
  outliers <- gps[which(gps$outlier == FALSE),]

  return(list(all = gps, retained = retained, outliers = outliers))
}
