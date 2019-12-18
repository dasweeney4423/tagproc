

load.argos <- function(argosfile, deployfile = NULL){
  argos <- read.csv(argosfile, stringsAsFactors = FALSE, row.names = NULL)
  if (names(argos)[1] == "row.names") {names(argos) <- c(names(argos)[2:ncol(argos)], "XNA")}
  if ("Satellite" %in% names(argos)) {argos <- argos[!argos$Satellite %in% c("13X0002", "14X0011"),]}
  if (!"LocationQuality" %in% names(argos)) {
    argos$LocationQuality <- as.character(argos$Quality)
    argos$LocationQuality[argos$Type %in% c("FastGPS", "User")] <- "3"
  }

  # Format time
  D <- time.turner(as.character(argos$Date))$strp
  argos$Date <- D

  # Check for deployment text file
  if (!is.null(deployfile)) {
    for (i in strsplit(argosfile, '/')[[1]]) {
      if (i == strsplit(argosfile, '/')[[1]][1]) {
        wcpath <- i
      } else {
        if (i == strsplit(argosfile, '/')[[1]][length(strsplit(argosfile, '/')[[1]])]) {
          break
        } else {
          wcpath <- paste(wcpath, i, sep = '/')
        }
      }
    }
    wclf <- list.files(wcpath)
    if (deployfile %in% wclf) {
      dep <- read.table(paste0(wcpath, "/", deployfile), sep=",", header=FALSE)
      if (nrow(dep) == 2) {dep <- dep[2,]}
      argos <- rbind(argos[1,], argos)
      ddt <- as.character(dep[1, 1])
      ddt <- time.turner(ddt)$strp
      argos$Date[1] <- ddt
      argos$Latitude[1] <- argos$Latitude2[1] <- as.numeric(as.character(dep[1, 2]))
      argos$Longitude[1] <- argos$Longitude2[1] <- as.numeric(as.character(dep[1, 3]))
      argos$LocationQuality[1] <- "3"
    }
  }

  return(argos)
}
