#' Plot spatial data on map
#'
#' This function is able to take many different spatial objects and plot them on an esri map.
#' @param ... Spatial objects that you wish to plot. Can be linestrings, points, polygons... Please make sure your projections are as you want them to be.
#' @param color A vector the same size as the number of spatial object inputs that specifies the color desired for each of the objects in the order they are given.
#' @param size A vector the same size as the number of spatial object inputs that specifies the size desired for each of the objects in the order they are given.
#' @param title A main title to be made for the figure
#' @param subtitle A subtitle to be made for the figure
#' @examples #examples not yet provided, sorry :(

spatial.map <- function(..., color = NULL, size = NULL, title = NULL, subtitle = NULL) {
  if (!is.null(color)) {
    if (length(color) != ...length()) {
      stop('number of color inputs must match the number of spatial object inputs')
    }
  }
  if (!is.null(size)) {
    if (length(size) != ...length()) {
      stop('number of size inputs must match the number of spatial object inputs')
    }
  }

  #get earth tiles
  esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

  #create all layers
  p <- ggplot2::ggplot() +
          ggspatial::annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none")

  for (l in 1:...length()) {
    data.layer <- list(...)[[l]]
    if ("SpatialPolygonsDataFrame" %in% class(data.layer)) {
      p <- p + ggspatial::layer_spatial(data = data.layer, size = size[l], color = color[l], fill = 'transparent')
    } else {
      p <- p + ggspatial::layer_spatial(data = data.layer, color = color[l], size = size[l])
    }
  }

  #include titles if desired
  if (!is.null(title)) {
    if (!is.null(subtitle)) {
      p <- p + ggplot2::ggtitle(title, subtitle = subtitle)
    } else {
      p <- p + ggplot2::ggtitle(title)
    }
  }

  #plot it all
  p
}
