time.since.obs <- function(env) {
  env$Date <- time.turner(as.character(env$Date))$strp
  env$locType <- as.character(env$locType)

  RESULT <- env
  NearestObs.Min <- rep(NA, times = nrow(env))

  obs <- env[env$locType == "o",]
  if (nrow(obs) > 0) {
    for (i in 1:nrow(env)) {
      # Calculate new location for StartTime
      time1 <- as.character(env$Date[i]) # Start time of ith  behavior record
      loci <- env$locType[i]
      if (loci == "o") {
        NearestObs.Min[i] <- 0
      } else {
        diffs <- difftime(obs$Date, time1, units = "mins") # Subtract all ARGOS times from time1
        priors <- which(diffs < 0) # Store all ARGOS times prior to time1 (negative)
        posts <- which(diffs > 0) # Store all ARGOS times after time1 (positive)
        if (length(priors) > 0 & length(posts) > 0) {
          prior <- priors[length(priors)] # Store most recent prior index
          post <- posts[1]  # Store earliest post index
          LastObs.Mins <- abs(round(as.numeric(diffs[prior]), digits = 2))
          NextObs.Mins <- round(as.numeric(diffs[post]), digits = 2)
          NearestObs.Min[i] <- min(c(LastObs.Mins, NextObs.Mins), na.rm = TRUE)
        }
      }
    }
  }
  RESULT$NearestObs.Min <- NearestObs.Min

  return(RESULT)
}
