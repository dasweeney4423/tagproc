prettify.reproc <- function(filepath, WCpretty_dir = NULL, DeployID,
                            dplat = NULL, dplong = NULL, dptime = NULL, dplc = 3, dperror = 100) {
  data <- read.csv(filepath, stringsAsFactors = FALSE, row.names = NULL)
  if (ncol(data) == 1) {
    data <- read.csv(filepath, stringsAsFactors = FALSE, row.names = NULL, sep = ';')
  }
  if (nrow(data) > 0) {
    colnames(data) <- c('Ptt', 'Latitude1', 'Longitude1', 'Latitude2', 'Longitude2',
                        'Quality',	'Date',	'Altitude',	'Pass',
                        'Satellite',	'Frequency',
                        'MsgCount', 'Error.radius',	'Error.Semi.major.axis', 'Error.Semi.minor.axis',
                        'Error.Ellipse.orientation','GDOP..m.Hz.','Dist..Subsat.Track..deg.',
                        'X')
    data$Longitude1 <- data$Longitude1 - 360
    data$Type <- 'Argos'

    if (!is.null(dplat)) {
      data <- rbind(NA, data)
      data$Ptt[1] <- data$Ptt[2]
      data$Latitude1[1] <- dplat
      data$Longitude1[1] <- dplong
      data$Date[1] <- dptime
      data$Quality[1] <- dplc
      data$Error.radius[1] <- dperror
      data$Type[1] <- 'User'
    }

    data$Date <- time.turner(data$Date)$raw
    names(data)[which(names(data) == "Quality")] <- "LocationQuality"
    data <- data[,c(1,2,3,6:ncol(data))]
    names(data)[which(names(data) == "Latitude1")] <- "Latitude"
    names(data)[which(names(data) == "Longitude1")] <- "Longitude"
    data$DeployID <- DeployID

    filename <- strsplit(filepath,'/')[[1]]
    filename <- filename[length(filename)]

    if (!is.null(WCpretty_dir)) {
      prettyname <- paste0(WCpretty_dir, "/", filename)
      write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
    } else {
      return(data)
    }
  }
}
