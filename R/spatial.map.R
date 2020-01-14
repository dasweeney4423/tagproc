spatial.map <- function(data, loc.color = 'black', loc.size = 1, title = NULL, subtitle = NULL,
                        spatial.size = NULL, spatial.color = NULL, spatial.fill = NULL, ...) {
  #Make sure the variable is projected as CRS = WGS84
  data.sp <- sp::SpatialPoints(data[c("Longitude", "Latitude")])
  proj4string(data.sp) = sp::CRS("+init=epsg:4326")
  data.sp <- sp::spTransform(data.sp, CRS("+init=epsg:4326"))

  #get earth tiles
  esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

  #create all layers
  suppressWarnings(suppressPackageStartupMessages(require(ggplot2)))
  p <- ggplot() +
          annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
          layer_spatial(data = zcdata_cat.sp, size = loc.size, color = loc.color)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (!is.null(subtitle)) {
    p <- p + ggtitle(subtitle = subtitle)
  }
  if (...length() > 0) {
    for (l in 1:...length()) {
      p <- p + layer_spatial(data = list(...)[[l]], size = spatial.size[l], color = spatial.color[l], fill = spatial.fill[l])
    }
  }

  #plot it all
  p
}
