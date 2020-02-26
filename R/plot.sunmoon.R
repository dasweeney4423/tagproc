#' Plot solar and lunar environmental data
#'
#' Plot environmental data that has been collected from the sun.moon() function
#' @param data A dataframe obtained after using the sun.moon() function that contains environmental data associated with tag data
#' @param xlab Optional input to specify desired x-axis label for the plot. Default is "Location Record"
#' @return A multipanelled figure plotting important environmental variables
#' @examples #examples not yet provided, sorry :(

plot.sunmoon <- function(data, xlab = "Location Record") {
  par(mfrow = c(4, 1))
  par(mar=c(.5,6,.5,.5), oma=c(4,2,2,2))
  cex.scale = 1

  plot(data$sunAzimuth, type = "o", pch = 16, cex = .5, ylim = c(0, 360), col = "orange", ann = FALSE, axes = FALSE)
  axis(2, las = 2, cex.axis = cex.scale)
  title(ylab = "Sun\nAzimuth", cex.lab = cex.scale)

  plot(data$moonAltitude, type = "o", pch = 16, cex = .5, col = "grey50", ann = FALSE, axes = FALSE)
  axis(2, las = 2, cex.axis = cex.scale)
  title(ylab = "Moon\nAltitude", cex.lab = cex.scale)

  plot(data$moonIlluminatedFraction, type = "o", pch = 16, cex = .5, col = "dodger blue", ann = FALSE, axes = FALSE)
  axis(2, las = 2, cex.axis = cex.scale)
  title(ylab = "Moon\nIlluminated\nFraction", cex.lab = cex.scale)

  plot(data$moonPhase, type = "o", pch = 16, cex = .5, col = "firebrick", ann = FALSE, axes = FALSE)
  axis(2, las = 2, cex.axis = cex.scale)
  axis(1, cex.axis = cex.scale)
  title(ylab = "Moon\nPhase", cex.lab = cex.scale)
  title(xlab = xlab, cex.lab = cex.scale, outer=T, line = 2)
}
