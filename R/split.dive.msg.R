#' Separate and process Behavior.csv data
#'
#' This functions takes data input from the Behavior.csv file from the WC portal and divides it into two dataframes, one for dive data and the other for data messages
#' @param data A dataframe containing the Behavior.csv data from the WC portal
#' @param kclusters The number of clusters desired to create for k-means clustering of dive depth and dive duration for the output of the Dives data. If left blank, no cluster will be performed.
#' @param lag The number of minutes for which sequence lags longer than it are considered data skips. Used when marking surfacings. Default is 15 minutes.
#' @return A list of two elements:
#' \itemize{
#' \item{\strong{Dives: }} A dataframe containing all surfacing or dive rows from the Behavior.csv data after having been processed
#' \item{\strong{Messages: }} A dataframe containing all message rows from the Behavior.csv data after having been processed
#' }
#' @examples #examples not yet provided, sorry :(

split.dive.msg <- function(data, kclusters = NULL, lag = 15) {
  # Format times (make sure seconds are kept)
  data$Start <- as.POSIXct(time.turner(as.character(data$Start))$strp)  # Start time
  data$End <- as.POSIXct(time.turner(as.character(data$End))$strp) # End time

  # Calculate duration of each row
  Duration <- difftime(time1 = data$End, time2 = data$Start, units = "mins")
  data$DurAvg <- round(as.numeric(Duration), digits = 4)

  # Separate messages and dives
  msg <- data[data$What == "Message",]
  dive <- data[data$What != "Message",]

  # Add record numbers
  msg <- msg[order(msg$Start),] # sort message data by start time
  msg$RecordNo <- c(1:nrow(msg)) # Add record number

  dive <- dive[order(dive$Start),] # sort message data by start time
  dive$RecordNo <- c(1:nrow(dive)) # Add record number

  # Calculate sequence lags
  # Dive data
  seqlag <- vector()
  seqlag[1] <- NA
  for (i in 2:nrow(dive)) {
    endi <- dive$End[i-1]
    starti <- dive$Start[i]
    lagi <- difftime(time1 = starti, time2 = endi, unit = "mins")
    seqlag[i] <- round(as.numeric(lagi), digits = 4)
  }
  dive$SeqLag <- seqlag

  # Message data
  if (nrow(msg) > 1) {
    seqlag <- vector()
    seqlag[1] <- NA
    for (i in 2:nrow(msg)) {
      endi <- msg$End[i-1]
      starti <- msg$Start[i]
      lagi <- difftime(time1 = starti, time2 = endi, unit = "mins")
      seqlag[i] <- round(as.numeric(lagi), digits = 4)
    }
  } else {
    seqlag <- 0
  }
  msg$SeqLag <- seqlag

  # Calculate average depth
  dive$DepthAvg <- round((dive$DepthMin + dive$DepthMax) / 2, digits = 4)

  #perform kmeans on dive data
  dive$Borderline <- dive$Kmeans <- NA
  surfacings <- dive[which(dive$What == 'Surface'),]
  if (!is.null(kclusters)) {
    kdata <- na.omit(data.frame(dive$DepthAvg, dive$DurAvg))
    kmdata <- na.omit(data.frame(scale(dive$DepthAvg), scale(dive$DurAvg)))
    k <- kmeans(kmdata, kclusters)
    kdata$kmeans <- k$cluster
    for (i in 1:nrow(kdata)) {
      diverow <- which(dive$DepthAvg == kdata$dive.DepthAvg[i] & dive$DurAvg == kdata$dive.DurAvg[i])
      dive$Kmeans[diverow] <- kdata$kmeans[i]
    }

    if (kclusters == 2) {
      if (mean(dive$DepthAvg[which(dive$Kmeans == 1)], na.rm=T) > mean(dive$DepthAvg[which(dive$Kmeans == 2)], na.rm=T)) {
        dive$Kmeans[which(dive$Kmeans != 1)] <- 'Shallow'
        dive$Kmeans[which(dive$Kmeans == 1)] <- 'Deep'
      } else {
        dive$Kmeans[which(dive$Kmeans != 1)] <- 'Deep'
        dive$Kmeans[which(dive$Kmeans == 1)] <- 'Shallow'
      }

      #mark the borderline clusterings
      deep <- dive[which(dive$Kmeans == 'Deep'),]
      deep$Borderline <- FALSE
      deep[which(scale(deep$DepthAvg) < quantile(scale(deep$DepthAvg), .05)),]$Borderline <- TRUE
      deep[which(scale(deep$DurAvg) < quantile(scale(deep$DurAvg), .05)),]$Borderline <- TRUE
      shallow <- dive[which(dive$Kmeans == 'Shallow'),]
      shallow$Borderline <- FALSE
      shallow[which(scale(shallow$DepthAvg) > quantile(scale(shallow$DepthAvg), .95)),]$Borderline <- TRUE
      shallow[which(scale(shallow$DurAvg) > quantile(scale(shallow$DurAvg), .95)),]$Borderline <- TRUE
      dive <- rbind(deep, shallow, surfacings)
      dive <- dive[order(dive$Start),]
    } else {
      dive$Borderline <- deep <- shallow <- NA
    }
  }

  # Format final output tables
  # Messages
  Messages <- data.frame(DeployID = msg$DeployID,
                         PTT = msg$Ptt,
                         RecordNo = msg$RecordNo,
                         Source = msg$Source,
                         Instr = msg$Instr,
                         Count = msg$Count,
                         Start = msg$Start,
                         End = msg$End,
                         SeqLag = msg$SeqLag,
                         What = msg$What,
                         DurationMin = msg$DurationMin,
                         DurationMax = msg$DurationMax,
                         DurAvg = msg$DurAvg)

  # Dives
  Dives <- data.frame(TagID = dive$DeployID,
                      PTT = dive$Ptt,
                      RecordNo = dive$RecordNo,
                      MsgCount = dive$Count,
                      StartTime = dive$Start,
                      EndTime = dive$End,
                      SeqLag = dive$SeqLag,
                      Event = dive$What,
                      Shape = dive$Shape,
                      DepthMin = dive$DepthMin,
                      DepthMax = dive$DepthMax,
                      DepthAvg = dive$DepthAvg,
                      DurationMin = dive$DurationMin,
                      DurationMax = dive$DurationMax,
                      DurAvg = dive$DurAvg,
                      Shallow = dive$Shallow,
                      Deep = dive$Deep,
                      Kmeans = dive$Kmeans,
                      Borderline = dive$Borderline)

  if (!is.null(kclusters)) {
    if (kclusters == 2) {
      Dives <- mark.surfacings(Dives, lag)
    }
  }

  return(list(Dives = Dives, Messages = Messages))
}
