#' Sync bathymetry and behavioral data
#'
#' Pulls bathymetry data from either the internet or a provided database and associates it with tag locations
#' @param data A dataframe of location data to associate with bathymetry data
#' @param z.radius Distance in km around each location that seafloor depth and slope will be calculated. This can be a single value or a vector of different values the same length as the number of rows of the data input
#' @param bathy.folder A directory path to the folder containing all bathymetry data. File names in the folder must be like this :'NEPACseafloor-116.000-30.000-115.000-31.000.csv'. If this input is NULL, data will be pulled from online, thus requiring internet access.
#' @return A dataframe with all tag data and its associated bathymetry data as new columns in the dataset
#' @examples #examples not yet provided, sorry :(

bathy.sync <- function(data, z.radius = 2.5, bathy.folder = NULL) {
  if (length(z.radius) > 1) {
    if (length(z.radius) != nrow(data)) {
      stop('length of z.radius must match number of rows in data or be a single, constant radius')
    }
  }
  if (!is.null(bathy.folder)) {
    # Prep data
    # Remove NA positions
    nas <- which(is.na(data[, which(names(data) %in% c("Latitude", "StartLat"))]))
    nas <- c(nas,which(is.na(data[, which(names(data) %in% c("Longitude", "StartLon"))])))
    nas <- unique(nas)
    if (length(nas) > 0) {data <- data[-nas,]}

    LAT <- data[, which(names(data) %in% c("Latitude", "StartLat"))]
    LON <- data[, which(names(data) %in% c("Longitude", "StartLon"))]

    latrange <- range(min(LAT, na.rm = TRUE), max(LAT, na.rm = TRUE))
    lonrange <- range(min(LON, na.rm = TRUE), max(LON, na.rm = TRUE))

    coord <- data.frame(longitude = LON, latitude = LAT)
    sp::coordinates(coord) <- c("longitude", "latitude")
    sp::proj4string(coord) <- sp::CRS("+proj=longlat +datum=WGS84")

    #########################################################
    # See which depth files will be needed for this tag
    lf <- list.files(bathy.folder)
    depthlf <- gsub("NEPACseafloor-", "", lf) # remove prefix
    depthlf <- gsub(".csv", "", depthlf) # remove ".csv"

    # Determine GPS boundaries of all the depth files
    lfsplit <- strsplit(depthlf, "-")
    lonmin <- lonmax <- latmin <- latmax <- vector()
    for (i in 1:length(lfsplit)) {
      lfspliti <- lfsplit[[i]]
      lonmin[i] <- as.numeric(paste0("-", lfspliti[1]))
      lonmax[i] <- as.numeric(paste0("-", lfspliti[3]))
      latmin[i] <- as.numeric(paste0(lfspliti[2]))
      latmax[i] <- as.numeric(paste0(lfspliti[4]))
    }

    # See which files include depth data for this tag
    toload <- vector()
    for (i in 1:length(depthlf)) {
      # Turn file into a polygon
      x1 <- lonmin[i]
      x2 <- lonmax[i]
      y1 <- latmin[i]
      y2 <- latmax[i]
      x <- c(x1, x1, x2, x2, x1)
      y <- c(y1, y2, y2, y1, y1)
      xy <- cbind(x, y)
      xy <- sp::Polygons(list(sp::Polygon(xy)), ID=1)
      xy <- sp::SpatialPolygons(list(xy), proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
      inpoly <- sp::over(coord, xy)
      nin <- length(inpoly[!is.na(inpoly)])
      if (nin > 0) {toload <- c(toload, i)}
    }

    # Load depth data
    zdf <- data.frame()
    for (i in 1:length(toload)) {
      zfile <- paste0(bathy.folder, '/', lf[toload[i]])
      zi <- read.csv(zfile, header = TRUE)
      zi <- zi[zi$z < 0,]
      zdf <- rbind(zdf, zi)
    }

    # Setup variables
    z <- zslope <- zaspect <- rep(NA, times = nrow(data))

    # Loop thru each entry, calculate variables
    for (i in 1:nrow(data)) {
      x <- LON[i]
      y <- LAT[i]
      if (length(z.radius) > 1) {
        z.rad <- z.radius[i]
      } else {
        z.rad <- z.radius
      }

      if (!is.na(x) & !is.na(y)) {
        # Seafloor depth
        top <- swfscMisc::destination(lat = y, lon = x, brng = 0, distance = z.rad, units = "km", type = "vincenty")[1]
        bottom <- swfscMisc::destination(lat = y, lon = x, brng = 180, distance = z.rad, units = "km", type = "vincenty")[1]
        left <- swfscMisc::destination(lat = y, lon = x, brng = 270, distance = z.rad, units = "km", type = "vincenty")[2]
        right <- swfscMisc::destination(lat = y, lon = x, brng = 90, distance = z.rad, units = "km", type = "vincenty")[2]
        zs <- zdf[zdf$x >= left & zdf$x <= right & zdf$y <= top & zdf$y >= bottom,]
        if (nrow(zs) > 0) {
          # Seafloor depth
          z[i] <- abs(round(mean(zs$z, na.rm = TRUE), digits = 3))

          # Seafloor slope
          zmaxi <- which.max(abs(zs$z))
          zmini <- which.min(abs(zs$z))
          zmax <- abs(zs$z[zmaxi])
          zmin <- abs(zs$z[zmini])
          zslop <- 100 * ((zmax - zmin) / (z.rad * 1000 * 2))
          zslope[i] <- round(zslop,digits = 3)

          # Get aspect of slope
          xmax <- zs$x[zmaxi]
          ymax <- zs$y[zmaxi]
          xmin <- zs$x[zmini]
          ymin <- zs$y[zmini]
          zaspect[i] <- as.numeric(swfscMisc::bearing(lat1 = ymax, lon1 = xmax, lat2 = ymin, lon2 = xmin)[1])
        }
      }
    }

    # Add calculated variables to results dataframe
    data$z <- z
    data$zslope <- zslope
    data$zaspect <- zaspect

    return(data)
  } else { #if online == TRUE
    # Prep data
    LAT <- data[, which(names(data) %in% c("Latitude", "StartLat"))]
    LON <- data[, which(names(data) %in% c("Longitude", "StartLon"))]
    if (z.radius < 2.5) {
      warning("This NOAA dataset can't use anything more precise than a radius of 2.5. If input is less than 2.5, it will be set to 2.5.")
    }

    # Get depth and slope for each point
    n <- zmin <- zmax <- z <- zslope <- zaspect <- rep(NA, times = nrow(data))
    for (i in 1:nrow(data)) {
      xi <- as.numeric(as.character(LON[i]))
      yi <- as.numeric(as.character(LAT[i]))
      if (length(z.radius) > 1) {
        z.rad <- z.radius[i]
      } else {
        z.rad <- z.radius
      }
      if (z.rad < 2.5) {
        z.rad <- 2.5 # This dataset can't use anything more precise than this radius
      }

      if (!is.na(xi) & !is.na(yi)) {
        top <- swfscMisc::destination(lat = yi, lon = xi, brng = 0, distance = z.rad, units = "km", type = "vincenty")[1]
        bottom <- swfscMisc::destination(lat = yi, lon = xi, brng = 180, distance = z.rad, units = "km", type = "vincenty")[1]
        left <- swfscMisc::destination(lat = yi, lon = xi, brng = 270, distance = z.rad, units = "km", type = "vincenty")[2]
        right <- swfscMisc::destination(lat = yi, lon = xi, brng = 90, distance = z.rad, units = "km", type = "vincenty")[2]
        bath <- suppressMessages(marmap::getNOAA.bathy(lon1 = left, lon2 = right, lat1 = bottom, lat2 = top, resolution = 1, keep = FALSE, antimeridian = FALSE))
        bath <- marmap::as.xyz(bath)
        names(bath) <- c("x", "y", "z")

        n[i] <- nrow(bath)
        if (min(bath$z, na.rm = T) >= 0) {
          zmin[i] <- min(bath$z, na.rm = T)
          zmax[i] <- max(bath$z, na.rm = T)
          z[i] <- mean(bath$z, na.rm = T)
        } else {
          zmin[i] <- abs(max(bath$z[bath$z < 0], na.rm = T))
          zmax[i] <- abs(min(bath$z, na.rm = T))
          z[i] <- abs(mean(bath$z[bath$z < 0], na.rm = T))
        }

        # Get slope - distance between points
        zdiff <- zmax[i] - zmin[i]
        lon1 <- bath$x[which.min(bath$z)]
        lon2 <- bath$x[which.max(bath$z)]
        lat1 <- bath$y[which.min(bath$z)]
        lat2 <- bath$y[which.max(bath$z)]
        zdist <- swfscMisc::distance(lat1, lon1, lat2, lon2, units = "km", method = "vincenty")
        zslope[i] <- zdiff / (zdist * 1000)
        zaspect[i] <- as.numeric(swfscMisc::bearing(lat1 = lat1, lon1 = lon1, lat2 = lat2, lon2 = lon2)[1])
      }
    }

    # Add calculated variables to results dataframe
    data$z <- z
    data$zslope <- zslope
    data$zaspect <- zaspect
    data$zn <- n
    data$zmin <- zmin
    data$zmax <- zmax

    return(data)
  }
}
