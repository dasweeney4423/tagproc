#' Find time cues for dives
#'
#' This function is used to find the time cues for the start and end of dives in a depth record and turn them into data with similar formatting to Wildlife Computers behavior log files.
#' @param tag A .nc file containing all datastreams from an archival tag or a tag data structure from load_nc()
#' @param mindepth The threshold in meters at which to recognize a dive or flight. Dives shallower than mindepth will be ignored.
#' @param mindur The threshold in seconds at which to recognize a dive. Dives that are not blow mindepth for at least mindur will be ignored.
#' @param divestart The threshold in meters at which a submergence might be considered a dive. If NULL, then wet/dry data is used for this threshold and the input for wetdry is required.
#' @param wetdry The wet/dry sensor threshold for determining when the tag is going dry or going wet. Default is 100.
#' @param kclusters The number of clusters desired to create for k-means clustering of dive depth and dive duration for the output of the Dives data. If left blank, no cluster will be performed.
#' @param acoustics The csv of acoustic audit data from the tag
#' @return A data frame with one row for each dive found. The data will have identical olumn names as can be found in a Ptt-Behavior.csv file from the WC portal.
#' @note This function utilizes the find_dives, load_nc, and fir_nodelay functions in the tagtools package (https://github.com/stacyderuiter/TagTools).
#' @export

simplify.archival2 <- function(tag, mindepth, mindur, divestart = NULL, wetdry = 100, kclusters = NULL, acoustics = NULL) {
  if (is.character(tag)) {
    #load nc file
    load_nc <- function(file, which_vars=NULL){
      if (!grepl('.nc', file)){
        file <- paste(file, '.nc', sep='')
      }
      file_conn <- ncdf4::nc_open(file)
      #get variable names present in this file
      vars <- names(file_conn$var)
      if (!is.null(which_vars)){
        vars <- vars[vars %in% which_vars]
      }
      #read in the variables one by one and store in a list
      X <- list()
      for (v in 1:length(vars)){
        #get metadata for variable v
        X[[v]] <- ncdf4::ncatt_get(file_conn, vars[v])
        # remove redundant name label
        X[[v]]$name <- NULL
        field_names <- names(X[[v]])
        # add the actual data matrix or vector
        X[[v]]$data <- ncdf4::ncvar_get(file_conn, vars[v])
        # make sure the sensor data is the first element of X[[v]]
        X[[v]] <- X[[v]][c('data', field_names)]
      }
      # entries of X should match variable names from netCDF file
      names(X) <- vars
      # get metadata and add it to X as "info"
      X$info <- ncdf4::ncatt_get( file_conn , 0 )
      ncdf4::nc_close(file_conn)
      class(X) <- c('animaltag', 'list')
      return(X)
    }
    tag <- load_nc(tag)
  } else {
    if (!is.list(tag)) {
      stop('tag should either be a tag data structure from load_nc or the file path as an input')
    }
  }

  #either use wet/dry if NULL or a depth value if given
  if (is.null(divestart)) {
    #get data from tag input
    wet <- tag$wet$data
    wfs <- tag$wet$sampling_rate
    if ('depth' %in% names(tag)) {
      p <- tag$depth$data
      pfs <- tag$depth$sampling_rate
    }
    if ('depth_corrected' %in% names(tag)) {
      p <- tag$depth_corrected$data
      pfs <- tag$depth_corrected$sampling_rate
    }
    if ('CorrectedDepth' %in% names(tag)) {
      p <- tag$CorrectedDepth$data
      pfs <- tag$CorrectedDepth$sampling_rate
    }

    #find cues (seconds) of becoming wet/dry and associate to find submergence cues
    if (length(unique(wet)) == 2) {
      wetdry <- 0.5
    }
    tdry <- which(c(0, diff(wet > wetdry)) > 0) / wfs
    delete <- c()
    for (i in 1:length(tdry)) {
      if (p[(pfs*tdry[i])] > 3) {
        delete <- c(delete, i)
      }
    }
    if (length(delete) > 0) {
      tdry <- tdry[-delete]
    }
    twet <- which(c(0, diff(wet > wetdry)) < 0) / wfs
    dryafter <- c()
    for (i in 1:length(twet)) {
      dryafter <- c(dryafter, tdry[which(tdry > twet[i])[1]])
    }
    ml <- min(c(length(dryafter), length(tdry)))
    subs <- data.frame(start = twet[1:ml], end = dryafter[1:ml])

    #remove dives that don't go deeper than mindepth for mindur
    subs$dur <- subs$end - subs$start
    subs <- subs[which(subs$dur >= mindur),]
    subs$max <- 0
    delete <- c()
    for (i in 1:nrow(subs)) {
      if (round(pfs*subs$end[i]) > length(p)) {
        delete <- c(delete, i)
        next
      }
      pdive <- p[round(pfs*subs$start[i]):round(pfs*subs$end[i])]
      if (max(pdive) < mindepth) {
        delete <- c(delete, i)
        next
      }
      subs$max[i] <- max(pdive)
      deepdur <- length(which(pdive >= mindepth)) / pfs
      if (deepdur < mindur) {
        delete <- c(delete, i)
      }
    }
    dives <- subs[-delete,]

    #determine dive shapes and convert to times
    dives$shape <- NA
    for (d in 1:nrow(dives)) {
      depth <- p[(pfs*dives$start[d]):(pfs*dives$end[d])]
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
      wd <- which(scale(deep$max) < quantile(scale(deep$max), .05))
      if (length(wd) > 0) {
        deep[wd,]$Borderline <- TRUE
      }
      wd <- which(scale(deep$dur) < quantile(scale(deep$dur), .05))
      if (length(wd) > 0) {
        deep[wd,]$Borderline <- TRUE
      }
      deep[wd,]$Borderline <- TRUE
      shallow <- dives[which(dives$Kmeans == 'Shallow'),]
      shallow$Borderline <- FALSE
      shallow[which(scale(shallow$max) > quantile(scale(shallow$max), .95)),]$Borderline <- TRUE
      shallow[which(scale(shallow$dur) > quantile(scale(shallow$dur), .95)),]$Borderline <- TRUE
      dives <- rbind(deep, shallow)
      dives <- dives[order(dives$start),]
    } else {
      dives$Borderline <- NA
    }
    dives$Shallow <- dives$Deep <- NA

    #gathering surfacing times
    surfacings <- data.frame()
    for (d in 1:(nrow(dives)-1)) {
      start <- dives$end[d]
      end <- dives$start[d+1]
      dur <- end - start
      max <- NA
      type <- 'Surface'

      # determine deep and shallow durations
      # shallow = dry (e.g. ziphius) or above certain depth (e.g. 2m for physalus)
      # deep = below shallow threshold but not deep AND long enough to count as dive (e.g. not 50m and 30s)
      #Shallow
      wd <- wet[(start*wfs):(end*wfs)]
      Shallow <- length(wd[which(wd > wetdry)]) / wfs #anything over wetdry is considered dry

      #Deep
      sdepth <- p[(start*pfs):(end*pfs)]
      Deep <- (length(sdepth < mindepth) / pfs) - Shallow

      shape <- Kmeans <- Borderline <- NA
      row <- data.frame(start, end, dur, max, shape, type, Borderline, Kmeans, Deep, Shallow)
      surfacings <- rbind(surfacings, row)
    }

    #combine dives and surfacings
    beh <- rbind(dives, surfacings)
    beh <- beh[order(beh$start),]

    #if surfacing happened before first dive, add it to data
    pbefore <- p[1:(pfs*beh$start[1])]
    pbeforedeepdur <- length(which(pbefore >= mindepth)) / pfs
    if ((tdry[1] < beh$start[1]) & (pbeforedeepdur >= mindur)) {
      start <- tdry[1]
      end <- beh$start[1]
      dur <- end - start
      max <- NA
      type <- 'Surface'

      # determine deep and shallow durations
      # shallow = dry (e.g. ziphius) or above certain depth (e.g. 2m for physalus)
      # deep = below shallow threshold but not deep AND long enough to count as dive (e.g. not 50m and 30s)
      #Shallow
      wd <- wet[(start*wfs):(end*wfs)]
      Shallow <- length(wd[which(wd > wetdry)]) / wfs #anything over wetdry is considered dry

      #Deep
      sdepth <- p[(start*pfs):(end*pfs)]
      Deep <- (length(sdepth < mindepth) / pfs) - Shallow

      shape <- Kmeans <- Borderline <- NA
      row <- data.frame(start, end, dur, max, shape, type, Borderline, Kmeans, Deep, Shallow)
      beh <- rbind(row, beh)
    }

    #if surfacing happened after last dive, add it to data
    pafter <- p[(pfs*beh$end[nrow(beh)]):length(p)]
    pafterdeepdur <- length(which(pafter >= mindepth)) / pfs
    tpwet <- twet[which(twet < (length(p) / pfs))]
    if ((tpwet[length(tpwet)] > beh$end[nrow(beh)]) & (pafterdeepdur >= mindur)) {
      start <- beh$end[nrow(beh)]
      end <- tpwet[length(tpwet)]
      dur <- end - start
      max <- NA
      type <- 'Surface'

      # determine deep and shallow durations
      # shallow = dry (e.g. ziphius) or above certain depth (e.g. 2m for physalus)
      # deep = below shallow threshold but not deep AND long enough to count as dive (e.g. not 50m and 30s)
      #Shallow
      wd <- wet[(start*wfs):(end*wfs)]
      Shallow <- length(wd[which(wd > wetdry)]) / wfs #anything over wetdry is considered dry

      #Deep
      sdepth <- p[(start*pfs):(end*pfs)]
      Deep <- (length(sdepth < mindepth) / pfs) - Shallow

      shape <- Kmeans <- Borderline <- NA
      row <- data.frame(start, end, dur, max, shape, type, Borderline, Kmeans, Deep, Shallow)
      beh <- rbind(beh, row)
    }
  } else {

    ###############################################################################
    ###############################################################################

    #use find_dives from tagtools package for dive times
    find_dives <- function(p, mindepth, sampling_rate = NULL, surface = 1, findall = 0) {
      searchlen <- 20         #how far to look in seconds to find actual surfacing
      dpthresh <- 0.25        #vertical velocity threshold for surfacing
      dp_lp <- 0.25           #low-pass filter frequency for vertical velocity
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
          if (!missing(findall) | ((!identical(ks0, empty)) & (!identical(ks1, empty)))) {
            k <- k + 1
            if (identical(ks0, empty)) {
              ton[k] <- 1
            } else {
              ton[k] <- max(tsurf[ks0])
            }
            if (identical(ks1, empty) ) {
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
      dp <- fir_nodelay(matrix(c(0, diff(p)), ncol = 1) * sampling_rate,
                        n, dp_lp / (sampling_rate / 2))
      #for each ton, look back to find last time whale was at the surface
      #for each toff, look forward to find next time whale is at the surface
      dmax <- matrix(0, length(ton), 2)
      for (k in 1:length(ton)) {
        ind <- ton[k] + (-round(searchlen * sampling_rate):0)
        ind <- ind[which(ind > 0)]
        ki = max(which(dp[ind] < dpthresh))
        if (length(ki) == 0 | is.infinite(ki)) {
          ki <- 1
        }
        ton[k] = ind[ki] ;
        ind <- toff[k] + (0:round(searchlen * sampling_rate))
        ind <- ind[which(ind <= length(p))]
        ki <- min(which(dp[ind] > -dpthresh))
        if (length(ki) == 0 | is.infinite(ki)) {
          ki <- 1
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

    fir_nodelay <- function(x, n, fc, qual='low', return_coefs = FALSE){
      # input checking
      # ================================================================
      # make sure x is a column vector or matrix
      if (!(sum(class(x) %in% c('matrix', 'vector')))){
        x <- as.matrix(x)
      }
      if (is.vector(x)) x <- as.matrix(x, nrow=length(x))

      # in case of multi-channel data, make sure matrix rows are samples and columns are channels
      if (dim(x)[2] > dim(x)[1]) x <- t(x)

      # make sure n is even to ensure an integer group delay
      n <- floor(n/2)*2


      # generate fir filter
      # ============================================================
      h <- signal::fir1(n=n,w=fc, type=qual)

      if (return_coefs){
        return(h)
      }else{ # carry out filtering

        # append fake samples to start and end of x to absorb filter delay
        # (output from these will be removed before returning result to user)
        nofsampling_rate <- floor(n/2)
        top_pad <- matrix(x[nofsampling_rate:2,], ncol=ncol(x))
        bot_pad <- matrix(x[(nrow(x)-1):(nrow(x)-nofsampling_rate),], ncol=ncol(x))
        x_pad <- rbind(top_pad, x, bot_pad)

        # filter the signal
        # ============================================================
        # apply filter to padded signal
        y <- apply(x_pad, MARGIN = 2, FUN=signal::filter, filt=h, nrow=nrow(x_pad))

        # account for filter offset (remove padding)
        y = y[n-1+(1:nrow(x)),]

        return(y)
      }
    }

    ###############################################################################
    ###############################################################################

    #get data from tag input
    if ('depth' %in% names(tag)) {
      p <- tag$depth$data
      pfs <- tag$depth$sampling_rate
    }
    if ('depth_corrected' %in% names(tag)) {
      p <- tag$depth_corrected$data
      pfs <- tag$depth_corrected$sampling_rate
    }
    if ('CorrectedDepth' %in% names(tag)) {
      p <- tag$CorrectedDepth$data
      pfs <- tag$CorrectedDepth$sampling_rate
    }

    #get dive times
    if (length(which(is.na(p))) > 0) {
      #work around missing data in Zica-20190111-173186
      data1 <- find_dives(p[1:(which(is.na(p))[1]-1)], mindepth, sampling_rate = pfs, surface = divestart, findall = 0)
      data2 <- find_dives(p[(which(is.na(p))[length(which(is.na(p)))]+1):length(p)], mindepth, sampling_rate = pfs, surface = divestart, findall = 0)
      data2$start <- data2$start + (which(is.na(p))[length(which(is.na(p)))]+1) - 1
      data2$end <- data2$end + (which(is.na(p))[length(which(is.na(p)))]+1) - 1
      fdout <- rbind(data1[c(1:(nrow(data1)-1)),], data2[c(2:nrow(data2)),])
    } else {
      fdout <- suppressWarnings(find_dives(p, mindepth, sampling_rate = pfs, surface = divestart, findall = 0))
    }
    if (p[(fdout$end[nrow(fdout)]*pfs)] > divestart) {
      fdout <- fdout[-nrow(fdout),]
    }
    fdout$dur <- fdout$end - fdout$start

    #filter out dives that aren't greater than mindepth for at least mindur
    delete <- c()
    for (i in 1:nrow(fdout)) {
      pdive <- p[(pfs*fdout$start[i]):(pfs*fdout$end[i])]
      deepdur <- length(which(pdive >= mindepth)) / pfs
      if (deepdur < mindur) {
        delete <- c(delete, i)
      }
    }
    if (!is.null(delete)) {
      dives <- fdout[-delete,]
    } else {
      dives <- fdout
    }

    #determine dive shapes and convert to times
    dives$shape <- NA
    for (d in 1:nrow(dives)) {
      depth <- p[(pfs*dives$start[d]):(pfs*dives$end[d])]
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
    dives$Shallow <- dives$Deep <- NA

    #gathering surfacing times
    surfacings <- data.frame()
    for (d in 1:(nrow(dives)-1)) {
      start <- dives$end[d]
      end <- dives$start[d+1]
      dur <- end - start
      max <- NA
      if (length(which(is.na(p))) > 0) {
        #work around missing data in Zica-20190111-173186
        if (length(unique(c(start:end, 474452:476527))) != length(c(start:end, 474452:476527))) {
          type <- 'MissingData'
          Shallow <- Deep <- shape <- Kmeans <- Borderline <- NA
          row <- data.frame(start, end, dur, max, shape, type, Borderline, Kmeans, Deep, Shallow)
          surfacings <- rbind(surfacings, row)
          next
        } else {
          type <- 'Surface'
        }
      } else {
        type <- 'Surface'
      }

      # determine deep and shallow durations
      # shallow = dry (e.g. ziphius) or above certain depth (e.g. 2m for physalus)
      # deep = below shallow threshold but not deep AND long enough to count as dive (e.g. not 50m and 30s)
      #Shallow
      pd <- p[(start*pfs):(end*pfs)]
      Shallow <- length(pd[which(pd < divestart)]) / pfs #anything above divestart is considered dry

      #Deep
      sdepth <- p[(start*pfs):(end*pfs)]
      Deep <- (length(sdepth < mindepth) / pfs) - Shallow

      shape <- Kmeans <- Borderline <- NA
      row <- data.frame(start, end, dur, max, shape, type, Borderline, Kmeans, Deep, Shallow)
      surfacings <- rbind(surfacings, row)
    }

    #combine dives and surfacings
    dives <- dives[,names(surfacings)]
    beh <- rbind(dives, surfacings)
    beh <- beh[order(beh$start),]

    #find dives with incomplete dives
    if (length(which(is.na(p))) > 0) {
      #work around missing data in Zica-20190111-173186
      data1 <- find_dives(p[1:(which(is.na(p))[1]-1)], mindepth, sampling_rate = pfs, surface = divestart, findall = 1)
      data2 <- find_dives(p[(which(is.na(p))[length(which(is.na(p)))]+1):length(p)], mindepth, sampling_rate = pfs, surface = divestart, findall = 1)
      data2$start <- data2$start + (which(is.na(p))[length(which(is.na(p)))]+1) - 1
      data2$end <- data2$end + (which(is.na(p))[length(which(is.na(p)))]+1) - 1
      dall <- rbind(data1[c(1:(nrow(data1)-1)),], data2[c(2:nrow(data2)),])
    } else {
      dall <- suppressWarnings(find_dives(p, mindepth, sampling_rate = pfs, surface = divestart, findall = 1))
    }

    #if surfacing happened before first dive, add it to data
    pbefore <- p[1:(pfs*beh$start[1])]
    pbeforedeepdur <- length(which(pbefore >= mindepth)) / pfs
    pminbefore <- which(pbefore < divestart)[1]
    if ((pminbefore < (pfs*beh$start[1])) & (pbeforedeepdur >= mindur)) {
      start <- dall$end[1]
      end <- beh$start[1]
      dur <- end - start
      max <- NA
      type <- 'Surface'

      # determine deep and shallow durations
      # shallow = dry (e.g. ziphius) or above certain depth (e.g. 2m for physalus)
      # deep = below shallow threshold but not deep AND long enough to count as dive (e.g. not 50m and 30s)
      #Shallow
      pd <- p[(start*pfs):(end*pfs)]
      Shallow <- length(pd[which(pd < divestart)]) / pfs #anything above divestart is considered dry

      #Deep
      sdepth <- p[(start*pfs):(end*pfs)]
      Deep <- (length(sdepth < mindepth) / pfs) - Shallow

      shape <- Kmeans <- Borderline <- NA
      row <- data.frame(start, end, dur, max, shape, type, Borderline, Kmeans, Deep, Shallow)
      beh <- rbind(row, beh)
    }

    #if surfacing happened after last dive, add it to data
    pafter <- p[(pfs*beh$end[nrow(beh)]):length(p)]
    pafterdeepdur <- length(which(pafter >= mindepth)) / pfs
    pminafter <- which(pafter < divestart)[length(which(pafter < divestart))] + length(p[1:(pfs*beh$end[nrow(beh)])])
    if ((pminafter > (pfs*beh$end[nrow(beh)])) & (pafterdeepdur >= mindur)) {
      start <- beh$end[nrow(beh)]
      end <- dall$end[nrow(dall)]
      dur <- end - start
      max <- NA
      type <- 'Surface'

      # determine deep and shallow durations
      # shallow = dry (e.g. ziphius) or above certain depth (e.g. 2m for physalus)
      # deep = below shallow threshold but not deep AND long enough to count as dive (e.g. not 50m and 30s)
      #Shallow
      pd <- p[(start*pfs):(end*pfs)]
      Shallow <- length(pd[which(pd < divestart)]) / pfs #anything above divestart is considered dry

      #Deep
      sdepth <- p[(start*pfs):(end*pfs)]
      Deep <- (length(sdepth < mindepth) / pfs) - Shallow

      shape <- Kmeans <- Borderline <- NA
      row <- data.frame(start, end, dur, max, shape, type, Borderline, Kmeans, Deep, Shallow)
      beh <- rbind(beh, row)
    }
  }

  #create return data
  TagID <- tag$info$depid
  if (length(strsplit(TagID, '-')[[1]]) == 3) {
    PTT <- strsplit(TagID, '-')[[1]][3]
  } else {
    PtT <- NA
  }
  MsgCount <- NA
  SeqLag <- 0
  RecordNo <- c(1:nrow(beh))
  StartSeconds <- beh$start
  EndSeconds <- beh$end
  StartTime <- beh$start + as.POSIXct(tag$info$dephist_device_datetime_start, format = "%d-%m-%Y %H:%M:%S", tz = 'UTC')
  EndTime <- beh$end + as.POSIXct(tag$info$dephist_device_datetime_start, format = "%d-%m-%Y %H:%M:%S", tz = 'UTC')
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
                       StartSeconds, EndSeconds,
                       StartTime, EndTime, SeqLag, Event,
                       Shape, DepthMin, DepthMax, DepthAvg,
                       DurationMin, DurationMax, DurAvg, Shallow,
                       Deep, Kmeans, Borderline)

  if (kclusters == 2) {
    output <- mark.surfacings(output)
  }

  #if SMRT with acoustics, determine which dives have foraging
  if (!is.null(acoustics)) {
    tagon <- as.POSIXct(tag$info$dephist_device_datetime_start, format = "%d-%m-%Y %H:%M:%S", tz = 'UTC')

    clicks <- read.csv(acoustics)
    clicks$EventStart <- clicks$EventStart %>% lubridate::ymd_hms(tz = "UTC")
    clicks$EventEnd <- clicks$EventEnd %>% lubridate::ymd_hms(tz = "UTC")
    clicks <- clicks[which(clicks$eventType == 'FD          '),]
    output$Foraging <- output$ClickEndTime <- output$ClickStartTime <- NA
    for (i in 1:nrow(dives)) {
      if (output$StartTime[i] > max(clicks$EventEnd)) {
        output$Foraging[i] <- "Unknown"
        next
      }
      if ((output$TagID[i] == "Zica-20191111-94810") & 
          (output$Start[i] >= 438636 & output$Start[i] <= 1033614)) {
        output$Foraging[i] <- "Unknown"
        next
      }
      diveclicks <- clicks[which(clicks$EventStart >= output$StartTime[i] &
                                   clicks$EventEnd <= output$EndTime[i]),]
      if (nrow(diveclicks) > 0) {
        output$Foraging[i] <- TRUE
        output$ClickStartTime[i] <- min(diveclicks$EventStart)
        output$ClickEndTime[i] <- max(diveclicks$EventEnd)
      } else {
        output$Foraging[i] <- FALSE
      }
    }
    output$ClickStartTime <- as.POSIXct(output$ClickStartTime, tz="UTC", origin="1970-01-01")
    output$ClickEndTime <- as.POSIXct(output$ClickEndTime, tz="UTC", origin="1970-01-01")
  }

  #put gps locations onto behaviors
  gps <- tag$GPS_position$data
  gps <- gps[which(tag$GPS_residual$data[,2] <= 35 & !is.na(gps[,3])),]
  output$LocationOffsetMins <- output$Longitude <- output$Latitude <- NA
  for (i in 1:nrow(output)) {
    dls <- abs(output$StartSeconds[i] - gps[,1])
    wdls <- which(dls == min(dls))[1]
    output$Latitude[i] <- gps[wdls,2]
    output$Longitude[i] <- gps[wdls,3]
    output$LocationOffsetMins[i] <- dls[wdls]/60
  }

  #clean up data where weird extra events popped in
  output <- output[which(output$DurAvg > 0),]
  output$RecordNo <- NA
  suppressPackageStartupMessages(suppressWarnings(require(tidyverse)))
  output <- output %>% distinct()
  output$RecordNo <- c(1:nrow(output))

  return(output)
}
