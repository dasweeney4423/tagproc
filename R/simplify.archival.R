#' Find time cues for dives
#'
#' This function is used to find the time cues for the start and end of dives in a depth record and turn them into data with similar formatting to Wildlife Computers behavior log files.
#' @param data A .nc file containing all datastreams from an archival tag
#' @param sampling_rate The sampling rate of the sensor data in Hz (samples per second).
#' @param mindepth The threshold in meters at which to recognize a dive or flight. Dives shallow or flights lower than mindepth will be ignored.
#' @param mindur The threshold in seconds at which to recognize a dive. Dives shorter than mindur will be ignored.
#' @param divestart The threshold in meters at which a tag might be considered. This is only used for Shallow and Deep durations.
#' @param surface (optional) The threshold in meters at which the animal is presumed to have reached the surface. Default value is 1. A smaller value can be used if the dive/altitude data are very accurate and you need to detect shallow dives/flights.
#' @param findall (optional) When 1 forces the algorithm to include incomplete dives at the start and end of the record. Default is 0 which only recognizes complete dives.
#' @param kclusters The number of clusters desired to create for k-means clustering of dive depth and dive duration for the output of the Dives data. If left blank, no cluster will be performed.
#' @param lag The number of minutes for which sequence lags longer than it are considered data skips. Used when marking surfacings. Default is 15 minutes.
#' @return A data frame with one row for each dive found. The data will have identical olumn names as can be found in a Ptt-Behavior.csv file from the WC portal.
#' @note This function utilizes the find_dives function in the tagtools package (https://github.com/stacyderuiter/TagTools).
#' @export

simplify.archival <- function(data, mindepth, mindur, divestart = NULL, surface = 1, findall = 0, kclusters = 2, lag = 15) {
  #find dives code from tagtools package
  find_dives <- function(p, mindepth, sampling_rate = NULL, surface = 1, findall = 0) {
    if (nargs() < 2) {
      stop("inputs for p and mindepth are required")
    }
    if (is.list(p)) {
      sampling_rate <- p$sampling_rate
      p <- p$data
      if (is.null(p)) {
        stop("p cannot be an empty vector")
      }
    } else {
      if (nrow(p) == 1) {
        p <- t(p)
      }
      if (is.null(sampling_rate)) {
        stop("sampling_rate is required when p is a vector")
      }
    }

    if (is.null(surface)) {
      surface <- 1          #maximum p value for a surfacing
    }

    if (is.null(findall)) {
      findall <- 0
    }

    searchlen <- 20         #how far to look in seconds to find actual surfacing
    dpthresh <- 0.25        #vertical velocity threshold for surfacing
    dp_lp <- 0.5           #low-pass filter frequency for vertical velocity
    #find threshold crossings and surface times
    tth <- which(diff(p > mindepth) > 0)
    tsurf <- which(p < surface)
    ton <- 0 * tth
    toff <- ton
    k <- 0
    empty <- integer(0)
    #sort through threshold crossings to find valid dive start and end points
    for (kth in 1:length(tth)) {
      if (all(tth[kth] > toff)) {
        ks0 <- which(tsurf < tth[kth])
        ks1 <- which(tsurf > tth[kth])
        if (findall || ((!identical(ks0, empty)) & (!identical(ks1, empty)))) {
          k <- k + 1
          if (identical(ks0, empty)) {
            ton[k] <- 1
          } else {
            ton[k] <- max(tsurf[ks0])
          }
          if (identical(ks1, empty)) {
            toff[k] <- length(p)
          } else {
            toff[k] <- min(tsurf[ks1])
          }
        }
      }
    }
    #truncate dive list to only dives with starts and stops in the record
    ton <- ton[1:k]
    toff <- toff[1:k]
    #filter vertical velocity to find actual surfacing moments
    n <- round(4 * sampling_rate / dp_lp)
    dp <- fir_nodelay(matrix(c(0, diff(p)), ncol = 1) * sampling_rate, n, dp_lp / (sampling_rate / 2))
    #for each ton, look back to find last time whale was at the surface
    #for each toff, look forward to find next time whale is at the surface
    dmax <- matrix(0, length(ton), 2)
    for (k in 1:length(ton)) {
      ind <- ton[k] + (-round(searchlen * sampling_rate):0)
      ind <- ind[which(ind > 0)]
      ki = which(dp[ind] < dpthresh)
      if (length(ki)==0) {
        ki <- 1
      }else{
        ki <- max(ki)
      }
      ton[k] = ind[ki] ;
      ind <- toff[k] + (0:round(searchlen * sampling_rate))
      ind <- ind[which(ind <= length(p))]
      ki <- which(dp[ind] > -dpthresh)
      if (length(ki)==0) {
        ki <- 1
      }else{
        ki = min(ki)
      }
      toff[k] <- ind[ki]
      dm <- max(p[ton[k]:toff[k]])
      km <- which.max(p[ton[k]:toff[k]])
      dmax[k, ] <- c(dm, ((ton[k] + km - 1) / sampling_rate))
    }
    #assemble output
    t0 <- cbind(ton,toff)
    t1 <- t0 / sampling_rate
    t2 <- dmax
    t <- cbind(t1, t2)
    t <- matrix(t[stats::complete.cases(t)], byrow = FALSE, ncol = 4)
    T <- data.frame(start = t[,1], end = t[,2],
                    max = t[,3], tmax = t[,4])
    return(T)
  }

  #gather data
  if (!is.list(data)) {
    stop('input for data must be a list with depth and info lists')
  }
  if ('depth' %in% names(data)) {
    p <- data$depth
  }
  if ('depth_corrected' %in% names(data)) {
    p <- data$depth_corrected
  }
  if ('CorrectedDepth' %in% names(data)) {
    p <- data$CorrectedDepth
  }
  tagon <- as.POSIXct(data$info$dephist_device_datetime_start, format = "%d-%m-%Y %H:%M:%S", tz = 'GMT')

  #get dive times
  dives <- find_dives(p, mindepth, sampling_rate = NULL, surface, findall)

  #filter by duration and convert to times
  dives$dur <- dives$end - dives$start
  dives <- dives[which(dives$dur >= mindur),]

  #determine dive shapes and convert to times
  dives$shape <- dives$Shallow <- dives$Deep <- NA
  for (d in 1:nrow(dives)) {
    depth <- p$data[(p$sampling_rate*dives$start[d]):(p$sampling_rate*dives$end[d])]
    divetime <- length(depth)
    bottomtime <- length(which(depth >= (max(depth) * .8)))
    if (bottomtime > (.5 * divetime)) {
      dives$shape[d] <- 'Square'
    } else {
      if (bottomtime <= (.2 * divetime)) {
        dives$shape[d] <- 'V'
      } else {
        dives$shape[d] <- 'U'
      }
    }
  }
  dives$type <- 'Dive'

  #perform kmeans clustering
  dives$Kmeans <- dives$Borderline <- NA
  kdata <- na.omit(data.frame(depth = scale(dives$max), dur = scale(dives$dur)))
  k <- kmeans(kdata, kclusters)
  dives$Kmeans <- k$cluster
  #label clusters if there are two
  if (kclusters == 2) {
    if (max(dives$max[which(dives$Kmeans == 1)]) > max(dives$max[which(dives$Kmeans == 2)])) {
      dives$Kmeans[which(dives$Kmeans == 1)] <- 'Deep'
      dives$Kmeans[which(dives$Kmeans == 2)] <- 'Shallow'
    } else {
      dives$Kmeans[which(dives$Kmeans == 2)] <- 'Deep'
      dives$Kmeans[which(dives$Kmeans == 1)] <- 'Shallow'
    }

    #mark the borderline clusterings
    deep <- dives[which(dives$Kmeans == 'Deep'),]
    deep$Borderline <- FALSE
    deep[which(scale(deep$max) < quantile(scale(deep$max), .05)),]$Borderline <- TRUE
    deep[which(scale(deep$dur) < quantile(scale(deep$dur), .05)),]$Borderline <- TRUE
    shallow <- dives[which(dives$Kmeans == 'Shallow'),]
    shallow$Borderline <- FALSE
    shallow[which(scale(shallow$max) > quantile(scale(shallow$max), .95)),]$Borderline <- TRUE
    shallow[which(scale(shallow$dur) > quantile(scale(shallow$dur), .95)),]$Borderline <- TRUE
    dives <- rbind(deep, shallow)
    dives <- dives[order(dives$start),]
  } else {
    dives$Borderline <- NA
  }

  #gathering surfacing times
  surfacings <- data.frame()
  for (d in 1:(nrow(dives)-1)) {
    start <- dives$end[d]
    end <- dives$start[d+1]
    dur <- end - start
    max <- tmax <- NA
    type <- 'Surface'

    # determine deep and shallow durations
    # shallow = dry (e.g. ziphius) or above certain depth (e.g. 2m for physalus)
    # deep = below shallow threshold but not deep or long enough to count as dive (e.g. not 50m and 30s)
    if (is.null(divestart)) {
      if ('wet' %in% names(data) == FALSE) {
        stop('wet/dry sensor required in data input if divestart is NULL')
      }

      #Shallow
      wd <- data$wet$data[(start * data$wet$sampling_rate):(end * data$wet$sampling_rate)]
      Shallow <- length(wd[which(wd > 50)]) / data$wet$sampling_rate #anything over 50 is considered dry

      #Deep
      Deep <- NA
    } else {
      if (divestart > mindepth) {
        stop('I think you mixed up divestart and mindepth. Divestart is for triggering when a dive might be happening and mindepth is the threshold that then determines that a dive is happening for sure.')
      }

      #Shallow
      wd <- p$data[(start * p$sampling_rate):(end * p$sampling_rate)]
      Shallow <- length(wd[which(wd < divestart)]) / p$sampling_rate

      #Deep
      Deep <- NA
    }

    shape <- Kmeans <- Borderline <- NA
    row <- data.frame(start, end, max, tmax, dur, Deep, Shallow, shape, type, Borderline, Kmeans)
    surfacings <- rbind(surfacings, row)
  }

  #combine dives and surfacings
  beh <- rbind(dives, surfacings)
  beh <- beh[order(beh$start),]
  beh$start <- tagon + beh$start
  beh$end <- tagon + beh$end

  #create return data
  TagID <- data$info$depid
  PTT <- MsgCount <- SeqLag <- NA
  RecordNo <- c(1:nrow(beh))
  StartTime <- beh$start
  EndTime <- beh$end
  Event <- beh$type
  Shape <- beh$shape
  DepthMin <- DepthMax <- DepthAvg <- beh$max
  DurationMin <- DurationMax <- beh$dur
  DurAvg <- beh$dur / 60
  Shallow <- beh$Shallow
  Deep <- beh$Deep
  Kmeans <- beh$Kmeans
  Borderline <- beh$Borderline
  output <- data.frame(TagID, PTT, RecordNo, MsgCount,
                       StartTime, EndTime, SeqLag, Event,
                       Shape, DepthMin, DepthMax, DepthAvg,
                       DurationMin, DurationMax, DurAvg, Shallow,
                       Deep, Kmeans, Borderline)

  if (kclusters == 2) {
    output <- mark.surfacings(output, lag)
  }

  return(output)
}
