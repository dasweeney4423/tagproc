#' Determine if land is between two points
#'
#' This function utilizes ptolemy package land data and determines if there is a body of land between two user-specified locations
#' @param loc1 Longitude and latitude (in that order) of first location
#' @param loc2 Longitude and latitude (in that order) of second location
#' @param land File path to the shapefile to be used for calculations as can be read by `sf::st_read(land)` or land spatial object from ptolemy package
#' @param ... Additional inputs to be passed to the function input for land. See the ptolemy package for possible inputs.
#' @return TRUE if land is between locations and FALSE otherwise
#' @examples #examples not yet provided, sorry :(

landbetween <- function(loc1, loc2, land, ...) {
  suppressPackageStartupMessages(require(sf))
  suppressPackageStartupMessages(require(ptolemy))

  #read in ptolemy land
  land <- suppressMessages(suppressWarnings(land(...)))

  #determine if land is between two points
  locs <- rbind(locs1, locs2)
  sf_locs <- st_as_sfc(paste0("LINESTRING(", locs[1,1], ' ', locs[1,2], ', ', locs[2,1], ' ', locs[2,2], ')')) %>%
    sf::st_set_crs(4326) %>%
    sf::st_transform(as.numeric(strsplit(st_crs(land)$input,':')[[1]][2]))
  ints <- st_intersection(sf_locs, land)
  if (length(ints) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
