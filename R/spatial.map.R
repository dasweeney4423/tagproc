spatial.map <- function(..., color = NULL, size = NULL, title = NULL, subtitle = NULL, crs = NULL) {

  #get earth tiles
  esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

  #create all layers
  p <- ggplot2::ggplot() +
          ggspatial::annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none")

  for (l in 1:...length()) {
    data.layer <- list(...)[[l]]
    if (!is.null(crs)) {
      data.layer <- sp::spTransform(data.layer, sp::CRS(crs))
    } else {
      if (l == 1) {
        crs1 <- sp::proj4string(data.layer)
      } else {
        if (sp::proj4string(data.layer) != crs1) {
          data.layer <- sp::spTransform(data.layer, sp::CRS(crs1))
        }
      }
    }

    if (class(data.layer) == "SpatialPoints") {
      p <- p + ggspatial::layer_spatial(data = data.layer, color = color[l], size = size[l])
    } else {
      if (class(data.layer) == "SpatialPolygonsDataFrame") {

        p <- p + ggspatial::layer_spatial(data = data.layer, size = size[l], color = color[l], spatial.fill = 'transparent')
      }
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
