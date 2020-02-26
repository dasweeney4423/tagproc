#' Plot bathymetry data
#'
#' Plot bathymetry data that has been collected from the bathy.sync()) function
#' @param data A dataframe obtained after using the bathy.sync() function that contains bathymetry data associated with tag data
#' @param xlab Optional input to specify desired x-axis label for the plot. Default is "Location Record"
#' @return A multipanelled figure plotting important bathymetry variables
#' @examples #examples not yet provided, sorry :(

plot.seafloor <- function(data, xlab = "Location Record") {
  par(mfrow=c(2, 1))

  # Seafloor depth
  par(mar = c(3, 6, .5, .5))
  plot(data$z, lwd = 3, ylim = c(max(data$z, na.rm = TRUE), 0), col = "steelblue1", type = "l", axes = FALSE, ann = FALSE)
  axis(1, cex.axis = 1)
  axis(2, las = 2, cex.axis = 1)
  title(ylab = "Depth (m)", cex.lab = 1, line = 4)

  # Seafloor slope
  par(mar = c(5, 6, .5, .5))
  plot(data$zslope, lwd = 3, ylim = c(0, max(data$zslope, na.rm = TRUE)), col = "purple", type = "l", axes = FALSE, ann = FALSE)
  axis(1, cex.axis = 1)
  axis(2, las = 2, cex.axis = 1)
  title(ylab = "% Slope", cex.lab = 1, line = 4)
  title(xlab = xlab, cex.lab = 1)
}
