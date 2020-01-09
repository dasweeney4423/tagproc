#' Plot environmental data
#'
#' Plot environmental data that has been collected from the sun.moon() function
#' @param data A dataframe obtained after using the sun.moon() function that contains environmental data associated with tag data
#' @return A multipanelled figure plotting important environmental variables
#' @examples #examples not yet provided, sorry :(

plot.sunmoon <- function(data) {
  par(mfrow = c(4, 1))
  par(mar = c(1, 6, .5, .5), oma = c(2, 2, 2, 2))
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
  title(ylab = "Moon\nPhase", xlab = "Dive record", cex.lab = cex.scale)
}