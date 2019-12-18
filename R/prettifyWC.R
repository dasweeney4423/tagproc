#' Clean up WC portal data
#'
#' This stage copies WC portal data to a new folder and formats column names and date/time information in the process. All subsequent stages of processing will draw on these “cleaned up” versions of the data. This way, the raw data will never be manipulated or replaced by TagProcApp functions.
#' @param infolder A file path to the directory where all the tag data CSV files are stored
#' @param outfolder A file path specifying the folder location and name where you would like all new, pretty files to be created. If NULL, a folder named "WC-processed" will be created within the infolder path.
#' @export
#' @examples #examples not yet provided, sorry :(

prettifyWC <- function(infolder, outfolder = NULL) {
  #create output folder name if it's input is NULL
  if (is.null(outfolder)) {
    outfolder <- paste(infolder, 'WC-processed', sep = '/')
  }

  # Determine WC-portal folder
  lf <- list.files(infolder)

  # Create WC-pretty folder
  dir.create(outfolder)

  # Create names of files you need to prettify
  filetags <- paste(c("All", "Argos", "Behavior", "Corrupt", "Histos",
                      "Locations", "MinMaxDepth", "RawArgos", "Series", "SeriesRange",
                      "Status", "Summary", "Argos-imade", "SST", "FastGPS"),
                    '.csv', sep = '')
  files2prettify <- c()
  PTT <- deployID <- NULL
  for (f in 1:length(lf)) {
    filef <- strsplit(lf[f], '-')[[1]]
    if (filef[length(filef)] %in% filetags) {
      if ((is.null(PTT)) & (length(filef) == 2)) {
        PTT <- filef[1]
      }
      if ((is.null(deployID)) & (length(filef) == 4)) {
        deployID <- filef[1]
      }
      files2prettify <- c(files2prettify, lf[f])
    }
  }

  ##########################################################
  # Scan for files in WC-portal folder OTHER than these files,
  # and copy them to the WC-pretty Files:
  copy.unchanged <- which(! lf %in% files2prettify)
  i <- 15
  for (i in copy.unchanged) {
    fromi <- paste0(infolder, "/", lf[i])
    toi <- paste0(outfolder, "/", lf[i])
    file.copy(from = fromi, to = toi, overwrite = TRUE)
  }

  ##########################################################
  # Now process each file you wish to prettify, if it exists
  for (FILE in files2prettify) {
    if (FILE == paste0(PTT, "-All.csv")){
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        if (names(data)[1] == "row.names") {
          names(data) <- c("DeployID", names(data)[3:ncol(data)])
          data <- data[,1:(ncol(data) - 1)]
        }
        names(data)[which(names(data) == "Platform.ID.No.")] <- "Ptt"
        names(data)[which(names(data) == "Loc..date")] <- "Date"
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # ARGOS
    if (FILE == paste0(PTT, "-Argos.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        if (names(data)[1] == "row.names") {
          names(data) <- c("DeployID", names(data)[3:ncol(data)])
          data <- data[,1:(ncol(data) - 1)]
        }
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # BEHAVIOR
    if (FILE == paste0(PTT, "-Behavior.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        data$Start <- time.turner(data$Start)$raw
        data$End <- time.turner(data$End)$raw
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # CORRUPT
    if (FILE == paste0(PTT, "-Corrupt.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0){
        # Format times
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file=prettyname, quote=FALSE, row.names=FALSE)
      }
      next
    }

    ##########################################################
    # HISTOS
    if (FILE == paste0(PTT, "-Histos.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # LOCATIONS
    if (FILE == paste0(PTT, "-Locations.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        data$Date <- time.turner(data$Date)$raw
        data$MsgCount <- data$Satellite <- NA
        names(data)[which(names(data) == "Quality")] <- "LocationQuality"
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # MINMAXDEPTH
    if (FILE == paste0(PTT, "-MinMaxDepth.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # RAWARGOS
    if (FILE == paste0(PTT, "-RawArgos.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        PassDate <- paste0(data$PassDate, " ", data$PassTime)
        PassDate <- time.turner(PassDate)$raw
        data$PassDate <- PassDate
        data$PassTime <- NULL
        MsgDate <- paste0(data$MsgDate, " ", data$MsgTime)
        MsgDate <- time.turner(MsgDate)$raw
        data$MsgDate <- MsgDate
        data$MsgTime <- NULL
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # SERIES
     if (FILE == paste0(PTT, "-Series.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        DATE <- paste0(data$Day, " ", data$Time)
        DATE <- time.turner(DATE)$raw
        data$Day <- DATE
        data$Time <- NULL
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # SERIESRANGE
    if (FILE == paste0(PTT, "-SeriesRange.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        data$Start <- time.turner(data$Start)$raw
        data$End <- time.turner(data$End)$raw
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data,file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # STATUS
    if (FILE == paste0(PTT, "-Status.csv")) {
      readname <- paste0(infolder,"/",FILE) ; readname
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        data$RTC <- time.turner(data$RTC)$raw
        data$Received <- time.turner(data$Received)$raw
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # SUMMARY
    if (FILE == paste0(PTT, "-Summary.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        if(names(data)[1] == "row.names") {
          names(data) <- c("DeployID" ,names(data)[3:ncol(data)])
          data <- data[,1:(ncol(data) - 1)]
        }
        # Format times
        data$EarliestXmitTime <- time.turner(data$EarliestXmitTime)$raw
        data$LatestXmitTime <- time.turner(data$LatestXmitTime)$raw
        data$EarliestDataTime <- time.turner(data$EarliestDataTime)$raw
        data$LatestDataTime <- time.turner(data$LatestDataTime)$raw
        prettyname <- paste0(outfolder, "/" , FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # ARGOS-imade
    if (FILE == paste0(PTT, "-Argos-imade.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # SST
    if (FILE == paste0(PTT, "-SST.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Format times
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(outfolder,"/",FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    # FastGPS
    if (FILE == paste0(PTT, "-FastGPS.csv")) {
      readname <- paste0(infolder, "/", FILE)
      data <- read.table(readname, header = TRUE, sep = ",",
                         skip = 3, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        # Change Name to DeployID
        names(data)[which(names(data) == "Name")] <- "DeployID"
        # Format times
        DATE <- paste0(data$Day, " ", data$Time)
        data$Day <- time.turner(DATE)$raw
        data$Time <- NULL
        prettyname <- paste0(outfolder, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
      next
    }

    ##########################################################
    if (!is.null(deployID)) {
      # 2/3-FastGPS
      if (FILE == paste0(deployID, "-", PTT, "-2-FastGPS.csv")) {
        readname <- paste0(infolder, "/", FILE)
        data <- read.table(readname, header = TRUE, sep = ",",
                           skip = 3, stringsAsFactors = FALSE, row.names = NULL)
        if (nrow(data) > 0) {
          # Change Name to DeployID
          names(data)[which(names(data) == "Name")] <- "DeployID"
          # Format times
          DATE <- paste0(data$Day, " ", data$Time)
          data$Day <- time.turner(DATE)$raw
          data$Time <- NULL
          data$InitTime <- time.turner(data$InitTime)$raw
          prettyname <- paste0(outfolder, "/", FILE)
          write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
        }
        next
      }
      ##########################################################
      # 2/3-Locations
      if (FILE == paste0(deployID, "-", PTT, "-2-Locations.csv")) {
        readname <- paste0(infolder, "/", FILE)
        data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
        if (nrow(data) > 0) {
          # Format times
          data$Date <- time.turner(data$Date)$raw
          prettyname <- paste0(outfolder, "/", FILE)
          write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
        }
        next
      }
    }
  } # End of FILE loop
}
