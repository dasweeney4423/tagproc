#' Calculate distance to land
#'
#' This tool takes location data and determines the distance between each location and the nearest land location from the given shapefile
#' @param data A dataframe containing the location data
#' @param land File path to the shapefile to be used for calculations as can be read by `sf::st_read(land)` or land spatial object from ptolemy package
#' @param crs Desired coordinate reference system. Input should match epsg code of land data
#' @param ... Additional inputs to be passed to the function input for land. See the ptolemy package for possible inputs.
#' @return A dataframe similar to the input for "data", but an extra column specifying the distance that location is to land (in km)
#' @examples #examples not yet provided, sorry :(

dist2land <- function(data, land, crs, ...) {
  suppressPackageStartupMessages(require(sf))
  suppressPackageStartupMessages(require(ptolemy))

  if (!is.character(land)) {
    #read in ptolemy land
    land <- suppressWarnings(land(...))
    allland <- data.frame()
    for (l in 1:length(land$geometry)) {
      allland <- rbind(allland, st_coordinates(land$geometry[[l]])[,c("X","Y")])
    }

    #convert tag locations to spatial coordinates
    data <- data %>% st_as_sf(coords = c('Longitude','Latitude')) %>%
      st_set_crs(as.numeric(strsplit(st_crs(land)$input,':')[[1]][2]))

    #calculate distances to land
    allland <- allland %>% st_as_sf(coords = c('X','Y'), crs = as.numeric(strsplit(st_crs(land)$input,':')[[1]][2]))
    proj4string(allland) <- CRS("+proj=longlat +init=epsg:3310")
    data$Dist2Land <- geosphere::dist2Line(p = st_coordinates(data$geometry),
                                           line = st_coordinates(allland)) / 1000
  } else {
    #read in land shape file
    land <- st_read(land)$geometry

    #convert tag locations to spatial coordinates
    data <- data %>% st_as_sf(coords = c('Longitude','Latitude')) %>%
      st_set_crs(st_crs(land))

    #calculate distances to land
    data$Dist2Land <- geosphere::dist2Line(p = st_coordinates(data$geometry),
                                           line = st_coordinates(land$geometry)[,c("X","Y")]) / 1000
  }

  return(data)
}
