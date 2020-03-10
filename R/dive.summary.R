#' Calculate dive summary statistics
#'
#' This function is used to calculate dive and surface behavior statistics and uses a forgaging depth threshold to label dive types (or uses kmeans if desired and provided in dive data)
#' @param dive Dive and surface behavior data dataframe obtained by the split.dive.msg() function.
#' @param summfile File path to where you might like the output data to be written as a .csv
#' @param forage.cutoff Indicates how to separate deep and shallow dives. Can either be a depth (in m) or can indicates 'kmeans' which will take kmeans analysis data from dive input to pull two classes of dives
#' @return A sigle row dataframe with all dive and surface statistics. When using this function to gather data across many individuals, simply use rbind() to add all rows together.
#' @examples #examples not yet provided, sorry :(

dive.summary <- function(dive, summfile = NULL, forage.cutoff = 600) {
  summ <- data.frame()
  dive$DepthAvg <- as.numeric(as.character(dive$DepthAvg))
  dive$Event <- as.character(dive$Event)

  ID <- dive$TagID[1]
  N.records <- nrow(dive)
  Tot.t <- round(as.numeric(difftime(dive$StartTime[nrow(dive)], dive$StartTime[1], units = "hours")), digits = 1)

  surf <- dive[dive$Event == "Surface",]
  dive <- dive[dive$Event == "Dive",]

  # Dive summaries
  N.dives <- nrow(dive[dive$Event == "Dive",])
  if ((forage.cutoff == 'kmeans') & ('Kmeans' %in% names(dive))) {
    deep <- which(dive$Kmeans == 'Deep')
    shallow <- which(dive$Kmeans == 'Shallow')
    N.deeps <- length(deep)
    N.shallow <- length(shallow)
  } else {
    if ((forage.cutoff == 'kmeans') & (('Kmeans' %in% names(dive) == FALSE) | (length(levels(dive$Kmeans)) > 2))) {
      stop('Kmeans is not a column in dive data or there are too many groupings (2 maximum)')
    } else {
      deep <- which(dive$DepthAvg >= forage.cutoff & dive$Event == "Dive")
      shallow <- which(dive$DepthAvg <= forage.cutoff & dive$Event == "Dive")
      N.deeps <- length(deep)
      N.shallow <- length(shallow)
    }
  }

  Start.date <- dive$StartTime[1]
  Stop.date <- dive$StartTime[nrow(dive)]

  Median.z.all <- round(median(dive$DepthAvg, na.rm = TRUE), digits = 1)
  Min.z.all <- round(min(dive$DepthAvg, na.rm = TRUE), digits = 3)
  Max.z.all <- round(max(dive$DepthAvg, na.rm = TRUE), digits = 3)

  Median.t.all <- round(median(dive$DurAvg, na.rm = TRUE), digits = 1)
  Min.t.all <- round(min(dive$DurAvg, na.rm = TRUE), digits = 3)
  Max.t.all <- round(max(dive$DurAvg, na.rm = TRUE), digits = 3)

  Median.z.deep <- round(median(dive$DepthAvg[deep], na.rm = TRUE), digits = 1)
  Min.z.deep <- round(min(dive$DepthAvg[deep], na.rm = TRUE), digits = 3)
  Max.z.deep <- round(max(dive$DepthAvg[shallow], na.rm = TRUE), digits = 3)

  Median.t.deep <- round(median(dive$DurAvg[deep], na.rm = TRUE), digits = 1)
  Min.t.deep <- round(min(dive$DurAvg[deep], na.rm = TRUE), digits = 3)
  Max.t.deep <- round(max(dive$DurAvg[deep], na.rm = TRUE), digits = 3)

  Median.z.shallow <- round(median(dive$DepthAvg[shallow], na.rm = TRUE), digits = 1)
  Min.z.shallow <- round(min(dive$DepthAvg[shallow], na.rm = TRUE), digits = 3)
  Max.z.shallow <- round(max(dive$DepthAvg[shallow], na.rm = TRUE), digits = 3)

  Median.t.shallow <- round(median(dive$DurAvg[shallow], na.rm = TRUE), digits = 1)
  Min.t.shallow <- round(min(dive$DurAvg[shallow], na.rm = TRUE), digits = 3)
  Max.t.shallow <- round(max(dive$DurAvg[shallow], na.rm = TRUE), digits = 3)

  # Surface time summaries
  Median.t.surf <- round(median(surf$DurAvg, na.rm = TRUE), digits = 1)
  Min.t.surf <- round(min(surf$DurAvg, na.rm = TRUE), digits = 3)
  Max.t.surf <- round(max(surf$DurAvg, na.rm = TRUE), digits = 3)
  N.surf <- nrow(surf)

  summ <- data.frame(ID, forage.cutoff, N.records, N.dives, N.deeps, N.shallow, N.surf,
                     Tot.t, Start.date, Stop.date,
                     Median.z.all, Min.z.all, Max.z.all,
                     Median.t.all, Min.t.all, Max.t.all,
                     Median.z.deep, Min.z.deep, Max.z.deep,
                     Median.t.deep, Min.t.deep, Max.t.deep,
                     Median.z.shallow, Min.z.shallow, Max.z.shallow,
                     Median.t.shallow, Min.t.shallow, Max.t.shallow,
                     Median.t.surf, Min.t.surf, Max.t.surf)

  if (!is.null(summfile)) {
    # Export tables
    write.csv(summ, file = summfile, quote = FALSE, row.names = FALSE)
  }

  return(summ)
}

