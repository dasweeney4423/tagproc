mark.surfacings <- function(data) {
  for (w in unique(data$DeployID)) {
    whale <- data[which(data$DeployID == w),]

  }
}
