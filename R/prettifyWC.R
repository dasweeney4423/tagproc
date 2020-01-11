#' Prettify files in the WC portal folder
#'
#' Copies WC portal data to a new folder and formats column names and date/time information in the process. All other data analyses will require these column headers.
#' @param WCfolder A file path to the WC portal folder with the tag data.
#' @return Once finished, a new folder names "WC-pretty" will be created within the folder specified by the input of 'WCfolder' in which all new data CSV files will be placed.
#' @export
#' @examples prettifyWC("WhaleTag001/WCPortal")

prettifyWC <- function(WCfolder) {
  #get files in WC-portal folder
  lf <- list.files(WCfolder)

  #determine PTT and deployID if present
  for (f in 1:length(lf)) {
    splitfile <- strsplit(strsplit(lf[f], '/')[[1]], '-')[[1]]
    if (splitfile[length(splitfile)] == 'All.csv') {
      PTT <- splitfile[1]
    }
    if (length(splitfile) > 1) {
      if ((splitfile[length(splitfile)] == 'FastGPS.csv') & (splitfile[(length(splitfile) - 1)] %in% c('2', '3'))) {
        deployID <- splitfile[1]
      }
    } else {
      deployID <- NULL
    }
  }

  #create WC-pretty folder
  newdir <- paste(WCfolder, 'WC-pretty', sep = '/')
  dir.create(newdir)

  #create names of files needing to be prettied
  filetags <- c("All", "Argos", "Behavior", "Corrupt", "Histos",
                "Locations", "MinMaxDepth", "RawArgos", "Series", "SeriesRange",
                "Status", "Summary", "Argos-imade", "SST", "FastGPS")
  files2prettify <- paste0(PTT, "-", filetags, ".csv")
  fastGPS <- paste0(deployID, "-", PTT, "-", c("2-FastGPS", "2-Locations", "3-FastGPS", "3-Locations"), ".csv")
  files2prettify <- c(files2prettify, fastGPS)
  files2prettify <- files2prettify[which(files2prettify %in% lf)]

  #process each file that needs to be prettied, if it exists
  for (FILE in files2prettify) {
    ##########################################################
    # ALL
    if (FILE == paste0(PTT, "-All.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        if (names(data)[1] == "row.names") {
          names(data) <- c("DeployID", names(data)[3:ncol(data)])
          data <- data[,1:(ncol(data) - 1)]
        }
        names(data)[which(names(data) == "Platform.ID.No.")] <- "Ptt"
        names(data)[which(names(data) == "Loc..date")] <- "Date"
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # ARGOS
    if (FILE == paste0(PTT, "-Argos.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        if (names(data)[1] == "row.names") {
          names(data) <- c("DeployID", names(data)[3:ncol(data)])
          data <- data[,1:(ncol(data) - 1)]
        }
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # BEHAVIOR
    if (FILE == paste0(PTT, "-Behavior.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        data$Start <- time.turner(data$Start)$raw
        data$End <- time.turner(data$End)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # CORRUPT
    if (FILE == paste0(PTT, "-Corrupt.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # HISTOS
    if (FILE == paste0(PTT, "-Histos.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # LOCATIONS
    if (FILE == paste0(PTT, "-Locations.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        data$Date <- time.turner(data$Date)$raw
        data$MsgCount <- data$Satellite <- NA
        names(data)[which(names(data) == "Quality")] <- "LocationQuality"
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # MINMAXDEPTH
    if (FILE == paste0(PTT, "-MinMaxDepth.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # RAWARGOS
    if (FILE == paste0(PTT, "-RawArgos.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data)>0) {
        PassDate <- paste0(data$PassDate, " ", data$PassTime)
        PassDate <- time.turner(PassDate)$raw
        data$PassDate <- PassDate
        data$PassTime <- NULL
        MsgDate <- paste0(data$MsgDate, " ", data$MsgTime)
        MsgDate <- time.turner(MsgDate)$raw
        data$MsgDate <- MsgDate
        data$MsgTime <- NULL
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names =FALSE)
      }
    }

    ##########################################################
    # SERIES
    if (FILE == paste0(PTT, "-Series.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if(nrow(data) > 0) {
        DATE <- paste0(data$Day, " ", data$Time)
        DATE <- time.turner(DATE)$raw
        data$Date <- DATE
        data$Time <- NULL
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }


    ##########################################################
    # SERIESRANGE
    if (FILE == paste0(PTT, "-SeriesRange.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        data$Start <- time.turner(data$Start)$raw
        data$End <- time.turner(data$End)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # STATUS
    if (FILE == paste0(PTT, "-Status.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname,stringsAsFactors=FALSE,row.names=NULL)
      if (nrow(data) > 0) {
        data$RTC <- time.turner(data$RTC)$raw
        data$Received <- time.turner(data$Received)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # SUMMARY
    if (FILE == paste0(PTT, "-Summary.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        if (names(data)[1] == "row.names") {
          names(data) <- c("DeployID", names(data)[3:ncol(data)])
          data <- data[,1:(ncol(data) - 1)]
        }
        data$EarliestXmitTime <- time.turner(data$EarliestXmitTime)$raw
        data$LatestXmitTime <- time.turner(data$LatestXmitTime)$raw
        data$EarliestDataTime <- time.turner(data$EarliestDataTime)$raw
        data$LatestDataTime <- time.turner(data$LatestDataTime)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # ARGOS-imade  # for FaRW031
    if(FILE==paste0(PTT,"-Argos-imade.csv")){
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # SST  # for GgTag018
    if (FILE == paste0(PTT, "-SST.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data)> 0) {
        data$Date <- time.turner(data$Date)$raw
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    ##########################################################
    # FastGPS  # for ZcTag053
    if (FILE == paste0(PTT, "-FastGPS.csv")) {
      readname <- paste0(WCfolder, "/", FILE)
      data <- read.table(readname, header = TRUE, sep=",",
                         skip = 3, stringsAsFactors = FALSE, row.names = NULL)
      if (nrow(data) > 0) {
        names(data)[which(names(data) == "Name")] <- "DeployID"
        DATE <- paste0(data$Day, " ", data$Time)
        data$Date <- time.turner(DATE)$raw
        data$Time <- NULL
        prettyname <- paste0(newdir, "/", FILE)
        write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
      }
    }

    if (!is.null(deployID)) {
      ##########################################################
      # 3/2-FastGPS  # for ZcTag053
      if (FILE %in% c(paste0(deployID, "-", PTT, "-2-FastGPS.csv"), paste0(deployID, "-", PTT, "-3-FastGPS.csv"))) {
        readname <- paste0(WCfolder, "/", FILE)
        data <- read.table(readname, header = TRUE, sep = ",",
                           skip = 3, stringsAsFactors = FALSE, row.names = NULL)
        if (nrow(data) > 0) {
          names(data)[which(names(data) == "Name")] <- "DeployID"
          DATE <- paste0(data$Day, " ", data$Time)
          data$Date <- time.turner(DATE)$raw
          data$Time <- NULL
          data$InitTime <- time.turner(data$InitTime)$raw
          prettyname <- paste0(newdir, "/", FILE)
          write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
        }
      }
      ##########################################################
      # 3/2-Locations  # for ZcTag053
      if (FILE %in% c(paste0(deployID, "-", PTT, "-2-Locations.csv"), paste0(deployID, "-", PTT, "-3-Locations.csv"))) {
        readname <- paste0(WCfolder, "/", FILE)
        data <- read.csv(readname, stringsAsFactors = FALSE, row.names = NULL)
        if (nrow(data) > 0) {
          data$Date <- time.turner(data$Date)$raw
          prettyname <- paste0(newdir, "/", FILE)
          write.csv(data, file = prettyname, quote = FALSE, row.names = FALSE)
        }
      }
    }
  }
}
