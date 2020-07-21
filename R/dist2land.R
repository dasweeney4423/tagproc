#' Calculate distance to land
#'
#' This tool takes location data and determines the distance between each location and the nearest land location from the given shapefile
#' @param data A dataframe containing the location data
#' @param land File path to the shapefile to be used for calculations as can be read by `sf::st_read(land)`
#' @return A dataframe similar to the input for "data", but an extra column specifying the distance that location is to land (in km)
#' @examples #examples not yet provided, sorry :(

dist2land <- function(data, land) {
  suppressPackageStartupMessages(require(sf))
  
  #read in land shape file
  land <- st_read(land)$geometry
  
  #convert tag locations to spatial coordinates
  data <- data %>% st_as_sf(coords = c('Longitude','Latitude')) %>% 
    st_set_crs(sf::st_crs(land))
  
  #calculate distances to land
  data$Dist2Land <- geosphere::dist2Line(p = st_coordinates(data$geometry), 
                                         line = st_coordinates(land)[,c("X","Y")]) / 1000
  
  return(data)
}