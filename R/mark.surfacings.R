#' Mark surfacing types
#'
#' This function is used to mark surfacings as either first, internal, terminal, only, or unknown
#' @param data A data frame as processed and returned by simplify.archival or split.dive.msg
#' @param lag The number of minutes for which sequence lags longer than it are considered data skips
#' @return A data frame similar to the data input with an extra column labeling the surface types
#' @export

mark.surfacings <- function(data, lag = 15) { #first, internal, terminal, unknown, only surfacing
  output <- data.frame()
  for (w in unique(data$TagID)) {
    whale <- data[which(data$TagID == w),]
    whale$SurfaceType <- NA
    deeps <- which(whale$Kmeans == 'Deep')

    #first surfacings
    firsts <- deeps + 1
    for (i in 1:nrow(whale)) {
      if (i %in% firsts == FALSE) {next}
      if (abs(whale$SeqLag[i]) > lag) {
        whale$SurfaceType[i] <- 'USS'
        next
      }
      if (i != nrow(whale)) {
        if ((abs(whale$SeqLag[i+1]) > lag)) {
          whale$SurfaceType[i] <- 'USS'
          next
        }
        if (whale$Kmeans[i+1] == 'Deep') {
          whale$SurfaceType[i] <- 'OSS'
        } else {
          whale$SurfaceType[i] <- 'FSS'
        }
      } else {
        whale$SurfaceType[i] <- 'USS'
      }
    }

    #terminal surfacings
    terminals <- deeps - 1
    for (i in 1:nrow(whale)) {
      if (i %in% terminals == FALSE) {next}
      if (!is.na(whale$SurfaceType[i])) {next}
      if (is.na(whale$SeqLag[i])) {
        whale$SurfaceType[i] <- 'USS'
        next
      }
      if (abs(whale$SeqLag[i]) > lag) {
        whale$SurfaceType[i] <- 'USS'
        next
      }
      if ((abs(whale$SeqLag[i+1]) > lag)) {
        whale$SurfaceType[i] <- 'USS'
        next
      }
      if (i < deeps[1]) {
        firstrows <- whale[which(whale$Event == 'Surface' & whale$RecordNo < deeps[1]),]
        if (nrow(firstrows) > 1) {
          whale$SurfaceType[i] <- 'TSS'
        } else {
          whale$SurfaceType[i] <- 'USS'
        }
      } else {
        whale$SurfaceType[i] <- 'TSS'
      }
    }

    #internal surfacings
    internals <- which(is.na(whale$SurfaceType) & whale$Event == 'Surface')
    first <- TRUE
    for (i in 1:nrow(whale)) {
      if (i %in% internals == FALSE) {next}
      if (!is.na(whale$SurfaceType[i])) {next}
      if (is.na(whale$SeqLag[i])) {
        whale$SurfaceType[i] <- 'USS'
        next
      }
      if (abs(whale$SeqLag[i]) > lag) {
        whale$SurfaceType[i] <- 'USS'
        next
      }
      if (i != nrow(whale)) {
        if ((abs(whale$SeqLag[i+1]) > lag)) {
          whale$SurfaceType[i] <- 'USS'
          next
        }
        if (i < deeps[1]) {
          firstrows <- whale[which(whale$Event == 'Surface' & whale$RecordNo < deeps[1]),]
          if (nrow(firstrows) >= 3) {
            if (first == FALSE) {
              whale$SurfaceType[i] <- 'ISS'
            } else {
              whale$SurfaceType[i] <- 'USS'
              first <- FALSE
            }
          } else {
            whale$SurfaceType[i] <- 'USS'
          }
        } else {
          whale$SurfaceType[i] <- 'ISS'
        }
      } else {
        whale$SurfaceType[i] <- 'USS'
      }
    }

    output <- rbind(output, whale)
  }
  return(output)
}
