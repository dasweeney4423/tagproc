plot.polygons <- function(data, polypath) {
  polylf <- list.files(polypath)
  polylf <- polylf[polylf != "zips"]
  par(mfrow = c(6, 2))
  for (j in 1:length(polylf)) {
    # Load polygon j
    polyname <- paste0(polypath, "/", polylf[j])
    polj <- rgdal::readOGR(dsn = polyname)
    plot(polj, col = adjustcolor("black", alpha.f = .2), main = gsub(".kml", "", polylf[j]))
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
  }
}
