#' Plot animal track and kml polygons
#'
#' This tool can be used to plot animal tracks over spatial polygons. Assuming the input location data has been previously passed through the polyprox function, locations within each respective polygon will be marked in red.
#' @param data A dataframe containing the location data from which animal tracks will be plotted. If the data contains the output of the polyprox function, then locations that are within each respective polygon will be marked.
#' @param polypath File path to either a folder containing many kml files (the folder cannot contain any other types of files) or the path to a single kml file
#' @return A plot for each polygon provided in the inputs with animal track data overlayed if it is near the polygon boundaries. All plots are centered on the polygons, so if the animal locations are too far away from a given polygon, they may not show up on the plot.
#' @examples #examples not yet provided, sorry :(

plot.polygons <- function(data, polypath) {
  polylf <- list.files(polypath)

  if (length(polylf) == 0) {
    polyname <- polypath
    # Load polygon
    polj <- suppressWarnings(rgdal::readOGR(dsn = polyname, verbose = F))
    polyname <- strsplit(polyname,'/')[[1]][length(strsplit(polyname,'/')[[1]])]
    sp::plot(polj, col = adjustcolor("black", alpha.f = .2), main = polyname)
    lines(x = data$Longitude, y = data$Latitude)

    polyname <- gsub("-", "", polyname)
    polyname <- gsub(" ", "", polyname)
    corrcol <- vector()
    for (J in 1:ncol(data)) {
      corrcol[J] <- length(grep(names(data)[J], polyname))
    }
    corrcol <- which(corrcol > 0)
    corrcol <- corrcol[length(corrcol)]
    inx <- data$Longitude[which(data[, corrcol] == 0)]
    iny <- data$Latitude[which(data[, corrcol] == 0)]
    points(inx, iny, col = "firebrick", pch = 16, cex = .75)
  } else {
    # Loop through each Polygon
    for (j in 1:length(polylf)) {
      # Load polygon j
      polyname <- paste0(polypath, "/", polylf[j])
      polj <- suppressWarnings(rgdal::readOGR(dsn = polyname, verbose = F))
      polyname <- strsplit(polyname,'/')[[1]][length(strsplit(polyname,'/')[[1]])]
      sp::plot(polj, col = adjustcolor("black", alpha.f = .2), main = gsub(".kml", "", polylf[j]))
      lines(x = data$Longitude, y = data$Latitude)

      corrcol <- vector()
      for (J in 1:ncol(data)) {
        corrcol[J] <- length(grep(names(data)[J], polyname))
      }
      corrcol <- which(corrcol > 0)
      corrcol <- corrcol[length(corrcol)]
      inx <- data$Longitude[which(data[, corrcol] == 0)]
      iny <- data$Latitude[which(data[, corrcol] == 0)]
      points(inx, iny, col = "firebrick", pch = 16, cex = .75)
    }
  }
}
