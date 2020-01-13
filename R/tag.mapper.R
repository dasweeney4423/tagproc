#' Basic figure creation of track
#'
#' This function creates a simple figure showing the animal track locations as points along with a visualized path
#' @param data A dataframe containing locations to be plotted
#' @param zoom Amount of zooming in figure creation. Larger values zoom OUT more while small numbers zoom IN. Default is 1
#' @param lon.offset  Amount of shifting to occur along the longitude when creating the figure.
#' @param lat.offset Amount of shifting to occur along the latitude when creating the figure. Default is no shifting
#' @param lwd Line width. Default is 0.5
#' @param pch Shape of the points. Default is 1
#' @param cex Size of the points. Default is 0.5
#' @param col Color of points and lines for track
#' @param mars Margin dimensions of the desired plot.
#' @param mapfile Desired path and name (not including .pdf) of the figure file to be created. Default is that no image file is created.
#' @return A figure showing all locations is plotted and if mapfile is not NULL, a pdf image file is created in the location specified by the input mapfile
#' @examples #examples not yet provided, sorry :(

tag.mapper <- function(data, zoom = 1, lon.offset = 0, lat.offset = 0,
                       lwd = .5, pch = 1, cex = .5, col = "black", mars = c(1, 1, 1, 1),
                       mapfile = NULL) {
  data <- data[!is.na(data$Longitude) & !is.na(data$Latitude),]

  # Determine map boundaries
  xrange <- range(data$Longitude, na.rm = TRUE)
  xdiff <- (abs(xrange[1]) - abs(xrange[2]))
  xctr <- mean(xrange)
  xrad <-  xdiff / 2
  xoff <- lon.offset * xdiff
  xmin <- (xctr + xoff) - (xrad * zoom)
  xmax <- (xctr + xoff) + (xrad * zoom)

  yrange <- range(data$Latitude, na.rm = TRUE)
  ydiff <- (abs(yrange[2]) - abs(yrange[1]))
  yctr <- mean(yrange)
  yrad <- ydiff / 2
  yoff <- lat.offset * ydiff
  ymin <- (yctr + yoff) - (yrad * zoom)
  ymax <- (yctr + yoff) + (yrad * zoom)

  if (!is.null(mapfile)) {
    pdf(mapfile, width = 10, height = 10)
  }
  par(mar = mars)
  maps::map('mapdata::worldHires', fill = TRUE, col = "grey90", xlim = c(xmin, xmax), ylim = c(ymin, ymax), border = "grey50")
  axis(1, at = round(xrange, digits = 2))
  axis(2, at = round(yrange, digits = 2))
  points(x = data$Longitude, y = data$Latitude, cex = cex, pch = pch, col = col)
  lines(x = data$Longitude, y = data$Latitude, lwd = lwd, col = col)
  if (!is.null(mapfile)) {dev.off()}
}
