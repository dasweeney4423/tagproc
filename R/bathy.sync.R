bathy.sync <- function(env, z.radius = .25, online = FALSE) {
  if (online == FALSE) {
    # Prep data
    # Remove NA positions
    nas <- which(is.na(env[, which(names(env) %in% c("Latitude", "StartLat"))]))
    nas <- c(nas,which(is.na(env[, which(names(env) %in% c("Longitude", "StartLon"))])))
    nas <- unique(nas)
    if (length(nas) > 0) {env <- env[-nas,]}

    LAT <- env[, which(names(env) %in% c("Latitude", "StartLat"))]
    LON <- env[, which(names(env) %in% c("Longitude", "StartLon"))]

    latrange <- range(min(LAT, na.rm = TRUE), max(LAT, na.rm = TRUE))
    lonrange <- range(min(LON, na.rm = TRUE), max(LON, na.rm = TRUE))

    coord <- data.frame(longitude = LON, latitude = LAT)
    coordinates(coord) <- c("longitude", "latitude")
    proj4string(coord) <- CRS("+proj=longlat +datum=WGS84")

    #########################################################
    # See which depth files will be needed for this tag
    lf <- list.files("./Resources-Seafloor")
    depthlf <- gsub("NEPACseafloor-", "", lf) # remove prefix
    depthlf <- gsub(".csv", "", depthlf) # remove ".csv"

    # Determine GPS boundaries of all the depth files
    lfsplit <- strsplit(depthlf, "-")
    lonmin <- lonmax <- latmin <- latmax <- vector()
    for (i in 1:length(lfsplit)) {
      lfspliti <- lfsplit[[i]]
      lonmin[i] <- as.numeric( paste0("-", lfspliti[1]))
      lonmax[i] <- as.numeric( paste0("-", lfspliti[3]))
      latmin[i] <- as.numeric( paste0(lfspliti[2]))
      latmax[i] <- as.numeric( paste0(lfspliti[4]))
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
      xy <- Polygons(list(Polygon(xy)), ID=1)
      xy <- SpatialPolygons(list(xy), proj4string = CRS("+proj=longlat +datum=WGS84"))
      inpoly <- over(coord, xy)
      nin <- length(inpoly[!is.na(inpoly)])
      if (nin > 0) {toload <- c(toload, i)}
    }

    # Load depth data
    zdf <- data.frame()
    for (i in 1:length(toload)) {
      zfile <- paste0("./Resources-Seafloor/", lf[toload[i]])
      zi <- read.csv(zfile, header = TRUE)
      zi <- zi[zi$z < 0,]
      zdf <- rbind(zdf, zi)
    }

    # Setup variables
    z <- zslope <- zaspect <- rep(NA, times = nrow(env))

    # Loop thru each entry, calculate variables
    for (i in 1:nrow(env)) {
      x <- LON[i]
      y <- LAT[i]

      if (!is.na(x) & !is.na(y)) {
        # Seafloor depth
        top <- swfscMisc::destination(lat = y, lon = x, brng = 0, distance = z.radius, units = "km", type = "vincenty")[1]
        bottom <- swfscMisc::destination(lat = y, lon = x, brng = 180, distance = z.radius, units = "km", type = "vincenty")[1]
        left <- swfscMisc::destination(lat = y, lon = x, brng = 270, distance = z.radius, units = "km", type = "vincenty")[2]
        right <- swfscMisc::destination(lat = y, lon = x, brng = 90, distance = z.radius, units = "km", type = "vincenty")[2]
        zs <- zdf[zdf$x >= left & zdf$x <= right & zdf$y <= top & zdf$y >= bottom,]
        if (nrow(zs) > 0) {
          # Seafloor depth
          z[i] <- abs(round(mean(zs$z, na.rm = TRUE), digits = 3))

          # Seafloor slope
          zmaxi <- which.max(abs(zs$z))
          zmini <- which.min(abs(zs$z))
          zmax <- abs(zs$z[zmaxi])
          zmin <- abs(zs$z[zmini])
          zslop <- 100 * ((zmax - zmin) / (z.radius * 1000 * 2))
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
    env$z <- z
    env$zslope <- zslope
    env$zaspect <- zaspect

    return(env)
  } else { #if online == TRUE
    if (z.radius < 2.5) {z.radius <- 2.5} # This dataset can't use anything more precise than this radius

    # Prep data
    LAT <- env[, which(names(env) %in% c("Latitude", "StartLat"))]
    LON <- env[, which(names(env) %in% c("Longitude", "StartLon"))]

    # Get depth and slope for each point
    n <- zmin <- zmax <- z <- zslope <- zaspect <- rep(NA, times = nrow(env))
    for (i in 1:nrow(env)) {
      xi <- as.numeric(as.character(LON[i]))
      yi <- as.numeric(as.character(LAT[i]))

      if (!is.na(xi) & !is.na(yi)) {
        top <- swfscMisc::destination(lat = yi, lon = xi, brng = 0, distance = z.radius, units = "km", type = "vincenty")[1]
        bottom <- swfscMisc::destination(lat = yi, lon = xi, brng = 180, distance = z.radius, units = "km", type = "vincenty")[1]
        left <- swfscMisc::destination(lat = yi, lon = xi, brng = 270, distance = z.radius, units = "km", type = "vincenty")[2]
        right <- swfscMisc::destination(lat = yi, lon = xi, brng = 90, distance = z.radius, units = "km", type = "vincenty")[2]
        bath <- suppressMessages(getNOAA.bathy(lon1 = left, lon2 = right, lat1 = bottom, lat2 = top, resolution = 1, keep = FALSE, antimeridian = FALSE))
        bath <- as.xyz(bath)
        names(bath) <- c("x", "y", "z")

        n[i] <- nrow(bath)
        zmin[i] <- abs(max(bath$z[bath$z < 0]))
        zmax[i] <- abs(min(bath$z))
        z[i] <- abs(mean(bath$z[bath$z < 0]))

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
    env$z <- z
    env$zslope <- zslope
    env$zaspect <- zaspect
    env$zn <- n
    env$zmin <- zmin
    env$zmax <- zmax

    return(env)
  }
}
