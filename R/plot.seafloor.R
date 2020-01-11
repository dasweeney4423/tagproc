plot.seafloor <- function(data) {
  par(mfrow=c(2, 1))

  # Seafloor depth
  par(mar = c(3, 6, .5, .5))
  plot(data$z, lwd = 3, ylim = c(max(data$z, na.rm = TRUE), 0), col = "steelblue1", type = "l", axes = FALSE, ann = FALSE)
  axis(1, cex.axis = 1.5)
  axis(2, las = 2, cex.axis = 1.5)
  title(ylab = "Depth (m)", cex.lab = 1.5, line = 4)

  # Seafloor slope
  par(mar = c(5, 6, .5, .5))
  plot(data$zslope, lwd = 3, ylim = c(0, max(data$zslope, na.rm = TRUE)), col = "purple", type = "l", axes = FALSE, ann = FALSE)
  axis(1, cex.axis = 1.5)
  axis(2, las = 2, cex.axis = 1.5)
  title(ylab = "% Slope", cex.lab = 1.5, line = 4)
  title(xlab = "Dive record", cex.lab = 1.5)
}
