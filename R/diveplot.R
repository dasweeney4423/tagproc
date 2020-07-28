#' Plot dive profile using dive shapes
#'
#' This function is used to visualize part or all of a dive profile by plotting dives based on their individual dive shapes
#' @param data A dataframe that contains the dive depths and shapes from split.dive.msg()
#' @param tz Time zone for plotting. If NULL, GMT is used.
#' @param zmin Minimum depth limit for plot
#' @param zmax Maximum depth limit on plot. If NULL, the maximum depth of all dives is used
#' @param Re.begin Which dive to begin at? Default is begin at the first dive
#' @param Rec.end Which dive to end at? Default is end at the last dive
#' @param cex.axis Axis font size
#' @param lwd Line width
#' @param col Line color
#' @param mars Figure margins
#' @param title Figure title if desired.
#' @param plot.records Whether or not to plot record numbers for dives
#' @param plot.w Width of figure to be made when saving as jpeg
#' @param plot.ratio Ratio of width to height when saving figure as jpeg
#' @param plotfile File path and name of where to save figure as jpeg. Default is that a figure is not saved
#' @return A dive profile is plotted and saved as a jpeg is desired
#' @examples #examples not yet provided, sorry :(

diveplot <- function(data,
                     tz = NULL,
                     zmin = 0, zmax = NULL,
                     Rec.begin = NULL, Rec.end = NULL,
                     cex.axis = 1,
                     lwd = 2,
                     col = "black",
                     mars = c(5, 5, 1, 2),
                     title = NULL,
                     plot.records = FALSE,
                     plot.w = 12, plot.ratio = 2.25,
                     plotfile = NULL) {
  if (is.null(Rec.begin)) {
    Rec.begin <- 1
  }
  if (is.null(Rec.end)) {
    Rec.end <- nrow(data)
  }
  data <- data[Rec.begin:Rec.end,]
  data$TagID <- as.character(data$TagID)
  data$PTT <- as.character(data$PTT)
  data$Event <- as.character(data$Event)
  data$Shape <- as.character(data$Shape)
  data$Shape[is.na(data$Shape)] <- "su"

  mr <- data.frame()
  for (i in 1:nrow(data)) {
    tag <- data$TagID[i]
    PTT <- data$PTT[i]
    ID <- data$RecordNo[i]
    ev <- data$Event[i]
    z <- data$DepthAvg[i]
    if (ev == "Surface") {z <- 0}
    shape <- data$Shape[i]
    subev <- paste0(ev, "-start")
    hr <- data$StartTime[i]
    startrow <- data.frame(tag, PTT, ID, hr, z, ev = subev, shape)
    hr <- data$EndTime[i]
    subev <- paste0(ev, "-end")
    endrow <- data.frame(tag, PTT, ID, hr, z, ev = subev, shape)
    mr <- rbind(mr, startrow, endrow)
  }

  # Time diff
  tdiff <- vector()
  tdiff[1] <- 0
  for (i in 2:nrow(mr)) {
    time1 <- mr$hr[i-1]
    time2 <- mr$hr[i]
    tdiffi <- difftime(time2, time1, units = "mins")
    tdiff[i] <- as.numeric(tdiffi)
  }
  mr$tdiff <- tdiff

  # Cumulative time
  mr$tsum <- cumsum(tdiff)

  # Shape dive profiles
  ids <- unique(mr$ID)
  mrplot <- data.frame()
  for (i in 1:length(ids)) {
    matches <- which(mr$ID == ids[i])
    mri <- mr[matches,]
    mr1 <- mri[1,]
    mr1$z <- 0
    mr2 <- mri[2,]
    mr2$z <- 0

    type <- as.character(mri$shape[1])

    time1 <- mri$tsum[1]
    time2 <- mri$tsum[2]

    timediff <- time2 - time1
    z <- mr$z[matches[1]]
    transittime <- z / 3
    transittime <- transittime / 60

    if (type == "V") {
      newtime <- time1 + (timediff / 2)
      mri$tsum[1:2] <- newtime
    }

    if (type == "Square") {
      tnew1 <- time1 + transittime
      tnew2 <- time2 - transittime
      mri$tsum[1] <- tnew1
      mri$tsum[2] <- tnew2
    }

    if (type == "U") {
      # The U has 12 nodes
      mru <- mri
      mru <- mru[1,] # take start of dive
      mru <- rbind(mru, mru, mru, mru, mru, mru, mru, mru, mru, mru, mru, mru)
      zvec <- c(.9, .95, .97, .99, .995, 1)
      zvec <- c(zvec, rev(zvec))
      zvec <- zvec * z
      tvec <- c(0, .03, .06, .09, .12, .15)
      tvec <- transittime + timediff * tvec
      tvec <- c(time1 + tvec, time2 - rev(tvec))
      mru$z <- zvec
      mru$tsum <- tvec
      mri <- mru
    }

    mri <- rbind(mr1, mri, mr2)
    mrplot <- rbind(mrplot, mri)
  }

  # Determine depth scale
  if (is.null(zmax)) {zmax <- max(mr$z, na.rm = TRUE)}

  # Prepare X avis labels
  # Hour nations
  if (!is.null(tz)) {
    mrplot$hr <- lubridate::with_tz(mrplot$hr, tzone = tz)
  }
  totdiff <- mrplot$tsum[nrow(mrplot)]
  times <- seq(0, totdiff, length = 5)
  labs <- strptime(mrplot$hr[1], format = "%Y-%m-%d %H:%M:%S", tz = "GMT") + times * 60
  labs <- paste0(substr(labs, 6, 10), " ", substr(labs, 12, 16))

  # Label locations
  markers <- vector()
  for (i in 1:length(times)) {
    diffs <- abs(mrplot$tsum - times[i])
    markers[i] <- which.min(diffs)
  }
  ats <- mrplot$tsum[markers]

  # Remove obviously erroneous records
  mrplot$tsum[mrplot$tdiff > 200 & substr(mrplot$ev, 1, 3) == "Dive"] <- NA
  mrplot$tsum[mrplot$tdiff > 20  & substr(mrplot$ev, 1, 3) == "Surface"] <- NA

  # Plot it!
  if (!is.null(plotfile)) {
    plot.h <- plot.w / plot.ratio
    if (!grepl('\\.jpeg', plotfile)) {
      plotfile <- paste(plotfile, '.jpeg', sep='')
    }
    jpeg(filename = plotfile, width = plot.w, height = plot.h, units = 'in', quality = 100, res = 300)
  }
  par(mar = mars)
  plot(x = mrplot$tsum, y = mrplot$z, type = "l", ylim = c(zmax, zmin), lwd = lwd, col = col, ann = FALSE, axes = FALSE)
  axis(1, at = ats, lab = labs, cex.axis = cex.axis)
  axis(2, cex.axis = cex.axis)
  if (plot.records) {
    axis(3, at = ats, lab = mrplot$ID[markers], cex.axis = .6)
  }
  if (!is.null(tz)) {
    title(xlab = "Local Time", ylab = "Depth (m)", cex.lab = cex.axis)
  } else {
    title(xlab = "Time (GMT)", ylab = "Depth (m)", cex.lab = cex.axis)
  }
  if (!is.null(title)) {
    title(main = title)
  }
  if (!is.null(plotfile)) {
    dev.off()
  }
}

