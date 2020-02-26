#' Create Kernel Density Estimates of home ranges
#'
#' This function creates Kernel Density Estimations of home ranges using the adehabitatHR package
#' @param data A dataframe containing one or more animal tracks
#' @param perc A vector of desired percentiles for home range estimation
#' @param kernel.method Indicates the desired home range estimation smoothing parameter method to use. Default is the ad hoc method ('href')
#' @param unout Indicates the units of the home range output. Default is squared kilometers
#' @param process.tz The time zone of the animal that is used in time conversions during the processing phase (if desired). If no additional processing is desired, leave this as NULL. See Dalla Rosa et al 2008 for home range data processing procedure and reasoning
#' @return A list containing each of the home range estimation levels
#' @examples #examples not yet provided, sorry :(

home.range <- function(data, perc = c(50, 95), kernel.method = 'href', unout = 'km2', process.tz = NULL) {
  if (!is.null(process.tz)) {
    #convert Dates and times into class POSIXct Datetime and convert lc to characters
    if ('POSIXct' %in% class(data$Date) == FALSE) {
      data$Date <- as.POSIXct(as.character(data$Date), format = '%m/%d/%Y %H:%M:%S', tz='UTC')
    }

    #process data
    data$localtime <- lubridate::with_tz(data$Date, tzone = process.tz)
    data$local_yday <- as.character(as.POSIXlt(data$localtime)$yday + 1)
    data$local_year <- as.character(as.POSIXlt(data$localtime)$year + 1900)
    data$local_month <- as.character(as.POSIXlt(data$localtime)$mon + 1)
    data$local_day <- as.character(as.POSIXlt(data$localtime)$mday)
    pdata <- data.frame()
    #average location per day
    for (w in unique(data$DeployID)) {
      whale <- data[which(data$DeployID == w),]
      for (d in unique(whale$local_yday)) {
        day <- whale[which(whale$local_yday == d),]
        mLatitude <- mean(day$Latitude)
        mLongitude <- mean(day$Longitude)
        Date <- as.POSIXct(paste(day$local_year, day$local_month, day$local_day, sep="-"), format = "%Y-%m-%d", tz = process.tz)
        pdata <- rbind(pdata, data.frame(DeployID = day$DeployID[1], Date = Date, Latitude = mLatitude, Longitude = mLongitude))
      }
    }
    ppdata <- data.frame()
    #add locations during duty cycled over days
    for (w in unique(pdata$DeployID)) {
      whale <- pdata[which(pdata$DeployID == w),]
      for (d in 1:nrow(whale)) {
        if (d == nrow(whale)) {
          ppdata <- rbind(ppdata, whale[d,])
          break
        }
        if (diff(c(whale$Date[d],whale$Date[(d+1)])) == 1) {
          ppdata <- rbind(ppdata, whale[d,])
        } else {
          p1 <- c(whale$Longitude[d], whale$Latitude[d])
          p2 <- c(whale$Longitude[(d+1)], whale$Latitude[(d+1)])
          dm <- geosphere::distVincentyEllipsoid(p1, p2) / 2 #half the distance between locations
          bear <- geosphere::bearing(p1, p2)
          nloc <- geosphere::destPoint(p1, bear, dm) #estimated location for middle day
          ppdata <- rbind(ppdata, whale[d,], data.frame(DeployID = whale$DeployID[1], Date = (whale$Date[d] + 86400), Latitude = nloc[2], Longitude = nloc[1]))
        }
      }
    }
  } else {
    ppdata <- data
  }

  #First, create a spatial points variable for lon and Latitude so that R recognizes these two columns as geographic coordinates:
  data.sp <- sp::SpatialPoints(ppdata[c("Longitude", "Latitude")])
  #Make sure the variable is projected as CRS = WGS84
  sp::proj4string(data.sp) = sp::CRS("+init=epsg:4326")
  data.sp <- sp::spTransform(data.sp, sp::CRS("+init=epsg:4326"))

  #Estimate home range KDE using the getveticshr function:
  data.kde <- adehabitatHR::kernelUD(data.sp, h = kernel.method)
  kde <- list()
  for (i in 1:length(perc)) {
    kde[[i]] <- adehabitatHR::getverticeshr(data.kde, percent = perc[i], unin = "m", unout = "km2")
  }

  return(hr = kde)
}
