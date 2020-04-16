#' Associate location data with behavioral data
#'
#' Combines location data and behavioral data in order to associate locations with all desired behaviors for both the start and end of the behaviors
#' @param data A dataframe containing the behavioral data to which locations will be associated
#' @param locs A dataframe containing the locations/animal track from which all locations will be obtained
#' @param matching Logical. If TRUE, then it is assumed each row in locations matches to the times of the rows in data. Default is FALSE.
#' @return A dataframe similar to the input for "data", but extra columns will be attached containing associated location and movement data for the animal.
#' @examples #examples not yet provided, sorry :(

bhvr.interpolate <- function(data, locs, matching = FALSE) {
  # Format bhvr logs
  data$StartTime <- time.turner(data$StartTime)$strp
  data$EndTime <- time.turner(data$EndTime)$strp

  if (matching) {
    if (nrow(data) != nrow(locs)) {
      stop("row numbers don't match between inputs so can't match locations directly")
    }
    locs$Date <- time.turner(locs$Date)$strp
    mr <- cbind(data, locs)
    ddate <- which(names(mr) == "Date")
    mr <- mr[,-ddate]
    return(mr)
  }

  # Format locs positions
  if (names(locs)[1] == "row.names") {names(locs) <- c(names(locs)[2:ncol(locs)], "XNA")}
  if ("Satellite" %in% names(locs)) {locs <- locs[!locs$Satellite %in% c("13X0002", "14X0011"),]}
  locs$Date <- time.turner(locs$Date)$strp
  if (!"LocationQuality" %in% names(locs)) {locs$LocationQuality <- NA}

  # Write function:
  # Inputs: a time and the Argos locations dataframe
  argos.locate <- function(time1, locs) {
    if (!"MsgCount" %in% names(locs)) {locs$MsgCount <- NA}

    # For ALL locations (observed and predicted)
    diffs <- difftime(locs$Date, time1, units = "mins") # Subtract all ARGOS times from time1
    priors <- which(diffs < 0) # Store all ARGOS times prior to time1 (negative)
    posts <- which(diffs > 0) # Store all ARGOS times after time1 (positive)

    lat <- lon <- NA
    LastPos.Mins <- NextPos.Mins <- NA
    LastPos.Type <- NextPos.Type <- NA
    LastPos.Msg <- NextPos.Msg <- NA
    LastPos.Satellite <- NextPos.Satellite <- NA
    LastPos.LocationQuality <- NextPos.LocationQuality <- NA
    LastPos.Km <- NextPos.Km <- NA

    if (length(priors) > 0 & length(posts) > 0) { # If the length of each is at least 1, then we have a chance to determine location.
      prior <- priors[length(priors)] # Store most recent prior index
      post <- posts[1]  # Store earliest post index

      # Store details of nearest Argos fixes
      LastPos.Mins <- abs(round(as.numeric(diffs[prior]), digits = 2))
      NextPos.Mins <- round(as.numeric(diffs[post]), digits = 2)
      LastPos.Msg <- locs$MsgCount[prior]
      NextPos.Msg <- locs$MsgCount[post]
      if ("Satellite" %in% names(locs)){
        LastPos.Satellite <- as.character(locs$Satellite[prior])
        NextPos.Satellite <- as.character(locs$Satellite[post])
        LastPos.LocationQuality <- as.character(locs$LocationQuality[prior])
        NextPos.LocationQuality <- as.character(locs$LocationQuality[post])
      }

      # Determine prior and post locations
      priorloc <- c(locs$Longitude[prior], locs$Latitude[prior])
      postloc <- c(locs$Longitude[post], locs$Latitude[post])

      # Calculate bearing and speed between prior and post Argos readings
      if (!any(is.na(c(priorloc, postloc)))) {
        Bearing <- swfscMisc::bearing(lat1 = priorloc[2], lat2 = postloc[2], lon1 = priorloc[1], lon2 = postloc[1])[1] # Determine difference between prior and time1
        DistanceM <- geosphere::distVincentyEllipsoid(priorloc, postloc)
        ElapsedTime <- as.numeric(difftime(locs$Date[post], locs$Date[prior], units = "secs"))
        Speed = DistanceM / ElapsedTime

        # Calculate new location for time1 based on ARGOS travel speed and time difference from latest ARGOS record
        EventDiff <- as.numeric(difftime(time1, locs$Date[prior], units = "secs"))
        Distance <- EventDiff * Speed
        p2 <- geosphere::destPoint(priorloc, Bearing, Distance)
        lat <- p2[2]
        lon <- p2[1]

        # Store details of nearest Argos fixes
        if (!is.na(lat) && !is.na(lon)) {
          LastPos.Km <- swfscMisc::distance(lat1 = priorloc[2], lon1 = priorloc[1], lat2 = lat, lon2 = lon, units = "km", method = "vincenty")
          NextPos.Km <- swfscMisc::distance(lat1 = postloc[2], lon1 = postloc[1], lat2 = lat, lon2 = lon, units = "km", method = "vincenty")
        }
      }
    }

    return(data.frame(Lat = lat, Lon = lon,
                      LastPos.Mins, LastPos.Km, LastPos.Msg, LastPos.Satellite, LastPos.LocationQuality,
                      NextPos.Mins, NextPos.Km, NextPos.Msg, NextPos.Satellite, NextPos.LocationQuality))
  }

  # Use ARGOS records to interpolate location of behavioral records
  RESULT <- data.frame(stringsAsFactors = FALSE)
  for (i in 1:nrow(data)) {
    # Calculate new location for StartTime
    time1 <- data$StartTime[i] # Start time of ith  behavior record
    starts <- argos.locate(time1, locs)
    names(starts) <- paste0("Start", names(starts))

    # Calculate new location for EndTime
    time1 <- data$EndTime[i]
    ends <- argos.locate(time1, locs)
    names(ends) <- paste0("End", names(ends))

    # Combine Start and End dataframes
    resulti <- cbind(starts, ends)
    RESULT <- rbind(RESULT, resulti)
  }

  # Add columns to dive record dataframe
  MR <- cbind(data, RESULT)

  return(MR)
}
