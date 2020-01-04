#' Perform the freitas filtering method
#'
#' This function filters location data obtained from Argos using the Freitas et al. (2008) algorithm.
#' @param argos A dataframe of argos locations to be sent through the Freitas filer algorithm
#' @param vmax Speed threshold for filter method, in m/s. Default is 2 m/s
#' @param ang Angles of the spikes to be removed by freitas filter. Default is c(15, 25). No spikes are removed if ang=-1
#' @param distlim Lengths of the above spikes (see input 'ang'), in kilometers. Default is c(2.5, 5).
#' @return A dataframe of filtered locations
#' @examples #examples not yet provided, sorry :(

freitas.filter <- function(argos, vmax = 2, ang = c(15, 25), distlim = c(2.5, 5)) {
  argos <- argos[!is.na(argos$Latitude) & !is.na(argos$Longitude),]
  distlim <- distlim * 1000
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
