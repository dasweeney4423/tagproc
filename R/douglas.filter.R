#' Perform Douglas Argos Filter
#'
#' This funciton is used to perform the Douglas Argos-Filter Algorithm on animal-track data provided from the Argos System for one individual from one tag.
#' @param argos An R dataframe of the argos locations. The dataframe must only include argos locations for one individual and one tag. To filter multiple individuals and tags, this function must be run separately for each.
#' @param argos_method A string specifying the data is from CLS' least-squared method (input == 'LS') or from the Kalman method (input == 'K').
#' @param method A string specifying the desired Douglas Argos-filtering method to be used. Options include 'MRD', 'DAR', or 'HYB'.
#' @param keep_lc The lc for which argos locations equal to or better than are retained according to the filtering method specified. If not specified, no argos locations are unconditionally kept because of their lc.
#' @param maxredun A required input for all filtering methods. This variable (units = km) is used to filter the Argos locations such that every retained location has a temporally neear-consecutive location that is spatially within the defined distance. See the following notes for more information on how this input is used differently in each filtering method.
#' @param duplrec A string specifying what to do in the cases where two argos locations have the same time stamp. If the input is 'offset' (default), then replicate times are offset by one second. If the input is 'filter', then replicates are marked as outliers and removed from the data before filtering.
#' @param keeplast Logical. An input only required if using the MRD or HYB filtering method. If FALSE (default), the last argos location is not unconditionally retained during filtering. If TRUE, the last argos location is unconditionally retained regardless of whether or not is passes the other filtering requirements.
#' @param skiploc Logical. An input only required if using the MRD or HYB filtering method. If FALSE (default), locations that are initially marked as an outlier are retested in the following round and have the chance of being retained if passing the second round of testing. If FALSE, locations that are initially marked as outliers are not retested and remain listed as an outlier.
#' @param minrate An input only required if using the DAR or HYB filtering method. minrate (units = km/hr) should reflect an upper bound of sustainable movement rate over a period of hours, including potential assistance by winds or currents.
#' @param ratecoef An input only required if using the DAR or HYB filtering method. A constant used to set the minimum allowed angle between three locations. Larger values for ratecoef will be less tolerant of acute turning angles; choices are typically between 15 and 25.
#' @param r_only Logical. An input only required if using the DAR or HYB filtering method. If TRUE, then the test to see if the angle between three locations is less than the minimum allowed angle is skipped. If FALSE (default), this angle test is included in the filtering process.
#' @param xmigrate An input only required if using the HYB filtering method. This input is used to determine if migration events occurred between two consecutive locations that both passed the MRD filter and that have at least one locations that passed the DAR filter in between them. A migration event is said to occur if the distance between the two MRD-retained locations exceeds this input (xmigrate) multiplied by maxredun.
#' @param xoverrun An input only required if using the HYB filtering method. If a migration event is determined, this input is used to test all intervening locations. The distance between the starting location of the migration and each intervening location must be less than the distance between the start location and the end location of the migration plus xoverrun times maxredun. If this test is TRUE, three more tests are performed to see if the DAR-retained intervening locations should be retained.
#' @param xdirect An input only required if using the HYB filtering method. Used to test the directional deviation of a location during a determined migration event. The heading from the starting location of the migration to a location along the migration path must be within +/- xdirect degrees of the heading of the vector from the starting migration location to the end location of the migration.
#' @param xangle An input only required if using the HYB filtering method. Used to test the angular deviation of a location during a determined migration event. The angle formed between the start, intervening, and ending  location of a migration event must exceed xangle degrees.
#' @param xpercent An input only required if using the HYB filtering method. Used to test the distance deviation of a location during a determined migration event. The total length between the start and intervening migration locations and the intervening and end migrations locations must not exceed the length of just the start and end migration locations by more than xpercent percentage.
#' @param test_0a An input only required if using the HYB filtering method. If the location quality for an intervening location during a migration event is A or better, test_0a number of tests (see most recent inputs) need to be passed for the location to be retained.
#' @param test_bz An input only required if using the HYB filtering method.If the location quality for an intervening location during a migration event is Z or B, test_bz number of tests (see most recent inputs) need to be passed for the location to be retained.
#' @param best_of_day Logical. If TRUE, an additional filtering process is carried out the only retains the "best location of the day" according to the method specified by the input rankmeth. Default is FALSE.
#' @param minoffh An input only required if best_of_day is TRUE and if you desire to filter according to PTT duty cycles rather than by GMT days. This input (units = hours) should be a little bit greater than the minimum OFF duration of the duty cycle and a little bit less than the maximum ON duration of the duty cycle. See the notes section for more information about how to properly set this input.
#' @param rankmeth An input only required if best_of_day is TRUE. The best location for the given cluster (GMT day or duty cycle depending on if minoffh is specified or not) is chosen by the following hierarchy of DIAG variables. For an input of 1 (default): lc, iqx, iqy, and nbmes. For an input of 2: lc, iqx, nbmes, and iqy. For an input of 3: lc and nbmes. As the cluster of locations passes through each DIAG variable test, the location with the highest values according to the hierachy is retained and all other locations in that cluster are filtered out.
#' @return A list of three elements in similar format to the data provided in the input for argos (note: if the input for duplrec was 'offset', times might be slightly different that the original data; all lat/lon primary and alternate locations were initially checked to find the shortest animal route through al combinations of primary and alternate locations and the determined superior locations are given the column header 'lon' or 'lat'):
#' \itemize{
#' \item{\strong{all: }} A dataframe containing all input argos locations with an additional column specifying whether the Douglas Argos-Filter marked the given location as an outlier or not.
#' \item{\strong{retained: }} A dataframe containing all retained input argos locations following the execution of the Douglas Argos-Filter algorithm.
#' \item{\strong{outliers: }} A dataframe containing all filtered-out input argos locations following the execution of the Douglas Argos-Filter algorithm.
#' }
#' @note For more information about the Douglas Argos-Filter Algorithm, see https://www.movebank.org/node/38, where you can find links to Douglas et al. 2012 (Methods in Ecology and Evolution) and the Douglas Argos-Filter Algorithm manual.
#' @export
#' @examples #examples not yet provided, sorry :(

douglas.filter <- function(argos, argos_method, method, keep_lc = NULL, maxredun = NULL, duplrec = 'offset', #inputs required by all filtering methods
                           keeplast = FALSE, skiploc = FALSE, #inputs for MRD filter
                           minrate = NULL, ratecoef = NULL, r_only = FALSE, #inputs for DAR filter
                           xmigrate = NULL, xoverrun = NULL, xdirect = NULL, xangle = NULL, xpercent = NULL, test_0a = NULL, test_bz = NULL, #inputs for Best Hybrid filter
                           best_of_day = FALSE, minoffh = NULL, rankmeth = 1 #inputs for Best of Day post-filtering
) {
  #input checks
  if ((method != 'HYB') & (method != 'DAR') & (method != 'MRD')) {
    stop('Invalid input for method. Valid options are MRD, DAR, or HYB')
  }
  if (argos_method %in% c('LS', 'K') == FALSE) {
    stop('Invalid input for argos_method. Valid options are LS or K')
  }
  if (!is.data.frame(argos)) {
    stop('Input for argos must be a dataframe')
  }
  if (keep_lc %in% c(1:3,'A','B','Z') == FALSE) {
    stop('input for keep_lc is not a real, valid lc')
  }

  #create necessary objects for later
  outliers <- data.frame()

  #convert dates and times into class POSIXct datetime and convert lc to characters
  if ('POSIXct' %in% class(argos$Date) == FALSE) {
    argos$Date <- as.POSIXct(as.character(argos$Date), format = '%m/%d/%Y %H:%M:%S', tz='UTC')
  }
  argos$Quality <- as.character(argos$Quality)

  #temporarily change lc letters to numbers to allow for easier comparison
  user <- data.frame()
  for (r in 1:nrow(argos)) {
    if (argos$Quality[r] %in% c('User', 'DP')) {
      user <- rbind(user, argos[r,]) #pull out user defined locations to ensure they are kept after filtering
      argos$Quality[r] <- 4
    }
    if (argos$Quality[r] == 'Z') {
      argos$Quality[r] <- -3
    }
    if (argos$Quality[r] == 'B') {
      argos$Quality[r] <- -2
    }
    if (argos$Quality[r] == 'A') {
      argos$Quality[r] <- -1
    }
  }
  argos$Quality <- as.numeric(argos$Quality)
  if (nrow(user) > 0) {
    user$outlier <- FALSE
  }

  #handle instances where there are duplicate times for multiple rows of locations
  if (duplrec == 'offset') {
    if (all(duplicated(argos$Date) == FALSE)) { #if all locations are not duplicates, move on
      done <- TRUE
    } else {
      done <- FALSE
    }

    #fix duplicated rows
    while (done == FALSE) {
      argos$Date[duplicated(argos$Date)] <- argos$Date[duplicated(argos$Date)] + 1 #add one second to all duplicated datetimes

      #check again to see if duplicated are still present
      if (all(duplicated(argos$Date) == FALSE)) { #if all locations are not duplicates, move on
        done <- TRUE
      } else {
        done <- FALSE
      }
    }
  } else {
    if (duplrec == 'filter') {
      argos$outlier <- FALSE
      argos$outlier[duplicated(argos$Date)] <- TRUE

      #remove outliers from argos data
      outliers <- rbind(outliers, argos[which(argos$outlier == TRUE),])
      argos <- argos[which(argos$outlier == FALSE),]
    } else {
      stop('Invalid input for duplrec. Valid options are offset or filter.')
    }
  }

  ###############################################################################
  #########################primary vs. alternate#################################
  ###############################################################################

  #find best of primary and alternate locations if data in least-squares argos locations
  if (argos_method == 'LS') {
    Zclass <- argos[which(argos$Quality == -3),]
    Zclass$Latitude <- NA
    Zclass$Longitude <- NA
    user_rows <- argos[which(argos$Quality == -4),]
    user_rows$Latitude <- user_rows$Latitude_p
    user_rows$Longitude <- user_rows$Longitude_p
    argos <- argos[which(argos$Quality != -3),]

    #create starting locations and distances/loc-strings
    l1 <- 1
    l2 <- l1 + 1
    while (l2 <= nrow(argos)) {
      #get locations
      loc1p <- c(argos$Longitude_p[l1], argos$Latitude_p[l1])
      loc1a <- c(argos$Longitude_a[l1], argos$Latitude_a[l1])
      loc2p <- c(argos$Longitude_p[l2], argos$Latitude_p[l2])
      loc2a <- c(argos$Longitude_a[l2], argos$Latitude_a[l2])

      #find four distances from combinations
      distpp <- geosphere::distVincentyEllipsoid(loc1p, loc2p)
      distpa <- geosphere::distVincentyEllipsoid(loc1p, loc2a)
      distaa <- geosphere::distVincentyEllipsoid(loc1a, loc2a)
      distap <- geosphere::distVincentyEllipsoid(loc1a, loc2p)

      #create loc-strings
      if (l1 == 1) { #if this is the first time through the loop create the first two elements of the loc-string
        if (distpp <= distap) {
          winstrp <- c('p', 'p')
          windistp <- distpp
        } else {
          winstrp <- c('a', 'p')
          windistp <- distap
        }
        if (distpa <= distaa) {
          winstra <- c('p', 'a')
          windista <- distpa
        } else {
          winstra <- c('a', 'a')
          windista <- distaa
        }
      } else {
        #find cumulative distances from the best tracks thus far
        cumpp <- windistp + distpp
        cumpa <- windistp + distpa
        cumap <- windista + distap
        cumaa <- windista + distaa

        if (cumpp <= cumap) {
          winstrp <- c(winstrp, 'p')
          windistp <- cumpp
        } else {
          winstrp <- c(winstra, 'p')
          windistp <- cumap
        }
        if (cumpa <= cumaa) {
          winstra <- c(winstrp, 'a')
          windista <- cumpa
        } else {
          winstra <- c(winstra, 'a')
          windista <- cumaa
        }
      }

      #go to next set of locs
      l1 <- l2
      l2 <- l1 + 1
    }

    #get shortest string
    if (windistp < windista) {
      winstring <- winstrp
    } else {
      winstring <- winstra
    }

    #go through argos dataframe and create lat and lon columns for the shortest path
    for (r in 1:nrow(argos)) {
      if (winstring[r] == 'p') {
        argos$Latitude[r] <- argos$Latitude_p[r]
        argos$Longitude[r] <- argos$Longitude_p[r]
      } else {
        argos$Latitude[r] <- argos$Latitude_a[r]
        argos$Longitude[r] <- argos$Longitude_a[r]
      }
    }

    #combine argos, user, and Zclass dataframes for next loop
    if ((nrow(user_rows) > 0) & (nrow(Zclass) > 0)) {
      argos <- rbind(argos, user_rows, Zclass)
      argos <- argos[order(argos$Date),]
    } else {
      if (nrow(user_rows) > 0) {
        argos <- rbind(argos, user_rows)
        argos <- argos[order(argos$Date),]
      } else {
        if (nrow(Zclass) > 0) {
          argos <- rbind(argos, Zclass)
          argos <- argos[order(argos$Date),]
        }
      }
    }

    #perform second pass with user and lc == Z rows included
    if (length(which(is.na(argos$Latitude))) > 0) {
      for (i in which(is.na(argos$Latitude))) {
        loc_p <- c(argos$Longitude_p[i], argos$Latitude_p[i])
        loc_a <- c(argos$Longitude_a[i], argos$Latitude_a[i])
        if (i == 1) {
          loc_after <- c(argos$Longitude[i+1], argos$Latitude[i+1])
          distp_after <- geosphere::distVincentyEllipsoid(loc_p, loc_after)
          dista_after <- geosphere::distVincentyEllipsoid(loc_a, loc_after)
          if (distp_after <= dista_after) {
            argos$Latitude[i] <- argos$Latitude_p[i]
            argos$Longitude[i] <- argos$Longitude_p[i]
          } else {
            argos$Latitude[i] <- argos$Latitude_a[i]
            argos$Longitude[i] <- argos$Longitude_a[i]
          }
        } else {
          if (i == nrow(argos)) {
            loc_before <- c(argos$Longitude[i-1], argos$Latitude[i-1])
            distp_before <- geosphere::distVincentyEllipsoid(loc_p, loc_before)
            dista_before <- geosphere::distVincentyEllipsoid(loc_a, loc_before)
            if (distp_before <= dista_before) {
              argos$Latitude[i] <- argos$Latitude_p[i]
              argos$Longitude[i] <- argos$Longitude_p[i]
            } else {
              argos$Latitude[i] <- argos$Latitude_a[i]
              argos$Longitude[i] <- argos$Longitude_a[i]
            }
          } else {
            loc_before <- c(argos$Longitude[i-1], argos$Latitude[i-1])
            distp_before <- geosphere::distVincentyEllipsoid(loc_p, loc_before)
            dista_before <- geosphere::distVincentyEllipsoid(loc_a, loc_before)
            loc_after <- c(argos$Longitude[i+1], argos$Latitude[i+1])
            distp_after <- geosphere::distVincentyEllipsoid(loc_p, loc_after)
            dista_after <- geosphere::distVincentyEllipsoid(loc_a, loc_after)
            dp <- distp_before + distp_after
            da <- dista_before + dista_after
            if (dp <= da) {
              argos$Latitude[i] <- argos$Latitude_p[i]
              argos$Longitude[i] <- argos$Longitude_p[i]
            } else {
              argos$Latitude[i] <- argos$Latitude_a[i]
              argos$Longitude[i] <- argos$Longitude_a[i]
            }
          }
        }
      }
    }
  }

  ###############################################################################
  #######################internal angle function#################################
  ###############################################################################

  int_ang <- function (locA, locB, locC) {
    bearBA <- geosphere::bearing(locB, locA)
    if (bearBA < 0) {
      bearBA <- 180 + (180 - abs(bearBA))
    }
    bearBC <- geosphere::bearing(locB, locC)
    if (bearBC < 0) {
      bearBC <- 180 + (180 - abs(bearBC))
    }
    alpha <- max(c(bearBA, bearBC)) - min(c(bearBA, bearBC))
    return(alpha)
  }

  ###############################################################################
  ###############################methods#########################################
  ###############################################################################

  #run MRD filtering method
  if ((method == 'MRD') | (method == 'HYB')) {
    MRD <- function(argos, maxredun, keep_lc, keeplast, skiploc) {
      if (is.null(maxredun)) {
        stop('maxredun is a required input.')
      }
      if (is.null(keep_lc)) {
        keep_lc <- 4 #if NULL, keep_lc is 4 so that only user given lcs are kept
      } else {
        if (is.character(keep_lc)) {
          if (keep_lc == 'Z') {
            keep_lc <- -3
          } else {
            if (keep_lc == 'B') {
              keep_lc <- -2
            } else {
              if (keep_lc == 'A') {
                keep_lc <- -1
              }
            }
          }
        }
      }

      #create indices to be used in loop
      rA <- 1
      rB <- rA + 1
      rC <- rB + 1

      #loop through argos data and filter out locations with distances greater than maxredun
      filtered <- rep(NA, nrow(argos))
      filt_rA <- FALSE
      while (rC <= nrow(argos)) {
        #find locations of three locations
        locA <- c(argos$Longitude[rA], argos$Latitude[rA])
        locB <- c(argos$Longitude[rB], argos$Latitude[rB])
        locC <- c(argos$Longitude[rC], argos$Latitude[rC])

        #calculate distances between three locations
        distAB <- geosphere::distVincentyEllipsoid(locA, locB) / 1000
        distBC <- geosphere::distVincentyEllipsoid(locB, locC) / 1000
        distAC <- geosphere::distVincentyEllipsoid(locA, locC) / 1000

        #mark the outlier location where necessary but keep locations that are good
        filtered[c(rA, rB, rC)] <- TRUE
        if (distAB < maxredun) {
          filtered[c(rA, rB)] <- FALSE
        }
        if (distBC < maxredun) {
          filtered[c(rB, rC)] <- FALSE
        }
        if (distAC < maxredun) {
          filtered[c(rA, rC)] <- FALSE
        }

        #go to next indices
        if ((skiploc == TRUE) & (filtered[rB] == TRUE)) {
          rA <- rC
          rB <- rA + 1
          rC <- rB + 1
          if ((rC > nrow(argos)) & is.na(filtered[nrow(argos)])) {
            rC <- nrow(argos)
            rB <- rC - 1
            rA <- rB - 1
            filt_rA <- TRUE
            filt_ra_index <- rA
          }
        } else {
          rA <- rB
          rB <- rA + 1
          rC <- rB + 1
        }
      }

      #in the case where we had to do one last iteration because we were about to skip the last argos location, ensure that the previosly-filtered-skipped location is kept as filtered
      if (filt_rA == TRUE) {
        filtered[filt_ra_index] <- TRUE
      }

      #apply filtering results to argos dataframe
      argos$outlier <- filtered

      #force locations with lc >= keep_lc to be retained
      argos$outlier[which(argos$Quality >= keep_lc)] <- FALSE

      #force last location to be retained if desired
      if (keeplast == TRUE) {
        argos$outlier[nrow(argos)] <- FALSE
      }

      return(argos)
    }
  }

  #run DAR filtering method
  if ((method == 'DAR') | (method == 'HYB')) {
    DAR <- function(argos, maxredun, ratecoef, minrate, keep_lc) {
      if (is.null(maxredun)) {
        stop('maxredun is a required input.')
      }
      if ((is.null(ratecoef)) | (is.null(minrate))) {
        stop('minrate and ratecoef are required inputs for the DAR method.')
      }
      if (is.null(keep_lc)) {
        keep_lc <- 4 #if NULL, keep_lc is 4 so that only user given lcs are kept
      } else {
        if (is.character(keep_lc)) {
          if (keep_lc == 'Z') {
            keep_lc <- -3
          } else {
            if (keep_lc == 'B') {
              keep_lc <- -2
            } else {
              if (keep_lc == 'A') {
                keep_lc <- -1
              }
            }
          }
        }
      }

      #create indices to be used in loop
      rA <- 1
      rB <- rA + 1
      rC <- rB + 1
      rD <- rC + 1

      #loop through argos data and filter out locations with distances greater than maxredun
      subargos <- argos
      argos$outlier <- NA
      iteration <- 1
      while (iteration <= 5) {
        subargos$outlier <- NA
        while (rD <= nrow(subargos)) {
          skip <- FALSE
          #carry out steps
          #step 1
          #find locations of A and B
          locA <- c(subargos$Longitude[rA], subargos$Latitude[rA])
          locB <- c(subargos$Longitude[rB], subargos$Latitude[rB])
          #calculate AB distance
          distAB <- geosphere::distVincentyEllipsoid(locA, locB) / 1000
          if (distAB < maxredun) {
            subargos$outlier[rB] <- FALSE
            skip <- TRUE
          }

          #step 2
          if (skip == FALSE) {
            if (subargos$Quality[rB] >= keep_lc) {
              subargos$outlier[rB] <- FALSE
              skip <- TRUE
            }
          }

          #step 3
          if (skip == FALSE) {
            if (r_only == FALSE) { #this step is only carried out if the input for r_only is FALSE
              #calculate alpha and minAlpha
              locC <- c(subargos$Longitude[rC], subargos$Latitude[rC])
              distBC <- geosphere::distVincentyEllipsoid(locB, locC) / 1000
              if ((distBC != 0) & (distAB != 0)) {
                alpha <- int_ang(locA, locB, locC) #see above for code
                minAlpha <- abs(-25 + ratecoef * log(min(distAB, distBC)))
                if (alpha < minAlpha) {
                  subargos$outlier[rB] <- TRUE
                  skip <- TRUE
                }
              }
            }
          }

          #step 4
          if (skip == FALSE) {
            #find rate from A to B
            hours <- (as.numeric(subargos$Date[rB]) - as.numeric(subargos$Date[rA])) / 60 / 60
            rateAB <- distAB / hours
            if (rateAB > minrate) {
              subargos$outlier[rB] <- TRUE
              skip <- TRUE
            }
          }

          #step 5
          if (skip == FALSE) {
            if (iteration == 5) {
              #find rate from B to C and calculate distBD and distCD
              hours <- (as.numeric(subargos$Date[rC]) - as.numeric(subargos$Date[rB])) / 60 / 60
              rateBC <- distBC / hours
              if (rateBC > minrate) {
                subargos$outlier[rB] <- TRUE
              }
            } else {
              #find rate from B to C and calculate distBD and distCD
              hours <- (as.numeric(subargos$Date[rC]) - as.numeric(subargos$Date[rB])) / 60 / 60
              rateBC <- distBC / hours
              locD <- c(subargos$Longitude[rD], subargos$Latitude[rD])
              distAC <- geosphere::distVincentyEllipsoid(locA, locC) / 1000
              distBD <- geosphere::distVincentyEllipsoid(locB, locD) / 1000
              distCD <- geosphere::distVincentyEllipsoid(locC, locD) / 1000
              if ((rateBC > minrate) & ((distAB + distBD) > (distAC + distCD))) {
                subargos$outlier[rB] <- TRUE
              }
            }
          }

          #move to next set of locations
          if (is.na(subargos$outlier[rB])) { #if location set did not trigger any of the steps, location B is retained
            subargos$outlier[rB] <- FALSE
          }
          if (subargos$outlier[rB] == TRUE) {
            rA <- rC
            rB <- rA + 1
            rC <- rB + 1
            rD <- rC + 1
          } else {
            rA <- rB
            rB <- rA + 1
            rC <- rB + 1
            rD <- rC + 1
          }
        }

        #apply filtering results to argos dataframe
        for (i in 1:nrow(subargos)) {
          if (is.na(subargos$outlier[i])) {
            subargos$outlier[i] <- FALSE
          }
          argos$outlier[which(argos$Date == subargos$Date[i])] <- subargos$outlier[i]
        }
        subargos <- subargos[which(subargos$outlier == FALSE),]

        #move to next iteration
        iteration <- iteration + 1
      }

      return(argos)
    }
  }

  if (method == 'HYB') {
    HYB <- function(argos, keep_lc, maxredun,
                    keeplast, skiploc, #inputs for MRD filter
                    minrate, ratecoef, r_only, #inputs for DAR filter
                    xmigrate, xoverrun, xdirect, xangle, xpercent, test_0a, test_bz) { #inputs for HYB filter

      #first run MRD filter and retain these locations unconditionally
      mrd_retained <- MRD(argos, maxredun, keep_lc, keeplast, skiploc)
      mrd_retained <- mrd_retained[which(mrd_retained$outlier == FALSE),]
      mrd_retained$method <- 'MRD'

      #then run DAR filter
      dar_retained <- DAR(argos, maxredun, ratecoef, minrate, keep_lc)
      dar_retained <- dar_retained[which(dar_retained$outlier == FALSE),]
      dar_retained$method <- 'DAR'

      #create full dataset
      ret <- rbind(dar_retained, mrd_retained)
      ret <- ret[order(ret$Date),]
      data <- data.frame()
      for (t in unique(ret$Date)) { #remove duplicate times that were saved by both filter methods
        dt <- ret[which(ret$Date == t),]
        if (nrow(dt) > 1) {
          data <- rbind(data, dt[which(dt$method == 'MRD'),])
        } else {
          data <- rbind(data, dt)
        }
      }

      #find start and end cues of possible migrations
      wmrd <- which(data$method=='MRD')
      dw <- diff(wmrd)
      migstart <- wmrd[which(dw > 1)] #data row at start of possible migration
      migend <- c()
      for (i in wmrd) { #data row at end of possible migration
        if (i %in% migrows) {
          migend <- c(migend, wmrd[which(wmrd == i) + 1])
        }
      }

      #determine if possible migrations are indeed migrations
      for (m in 1:length(migstart)) {
        s <- migstart[m]
        e <- migend[m]
        locstart <- c(data$Longitude[s], data$Latitude[s])
        locend <- c(data$Longitude[e], data$Latitude[e])
        migdist <- geosphere::distVincentyEllipsoid(locstart, locend) / 1000 #distance from start to end of possible migration
        if (migdist > (xmigrate * maxredun)) { #if TRUE, migration event occurred and now test all DAR lcoations in between
          #distance between Xo and Xn-1 must be less than the distance between Xo and Xn + xoverrun*maxredun
          locDARend <- c(data$Longitude[e-1], data$Latitude[e-1])
          distDAR <- geosphere::distVincentyEllipsoid(locstart, locDARend) / 1000 #distance between Xo and Xn-1
          if (distDAR < (migdist + xoverrun * maxredun)) { #if TRUE, three more tests are performed
            data$ndev <- 0
            for (i in (s + 1):(e - 1)) {
              loci <- c(data$Longitude[i], data$Latitude[i])
              #test 1: directional deviation
              beari <- geosphere::bearing(locstart, loci)
              bearse <- geosphere::bearing(locstart, locend)
              #the heading of the vector Xo-Xi must be within ±xdirect degrees of the heading of the vector Xo-Xn
              if (abs(beari - bearse) <= xdirect) {
                data$ndev[i] <- data$ndev[i] + 1
              }

              #test 2: angular deviation
              devang <- int_ang(locstart, loci, locend) #see above for code
              #the angle formed by Xo-Xi-Xn must exceed XANGLE degrees
              if (devang > xangle) {
                data$ndev[i] <- data$ndev[i] + 1
              }

              #test 3: distance deviation
              distsi <- geosphere::distVincentyEllipsoid(locstart, loci) / 1000
              distie <- geosphere::distVincentyEllipsoid(loci, locend) / 1000
              tdist <- distsi + distie
              #the length of Xo-Xi-Xn must not exceed the length of Xo-Xn by more than xpercent
              if (tdist <= (migdist + (migdist * xpercent))) {
                data$ndev[i] <- data$ndev[i] + 1
              }

              #retain or filter based on lc and other inputs
              if (data$Quality[i] %in% c(-3, -2)) {
                if (data$ndev[i] >= test_bz) {
                  data$outlier[i] <- FALSE
                }
              } else {
                if (data$ndev[i] >= test_0a) {
                  data$outlier[i] <- FALSE
                }
              }
            }
          }
        }
      }

      argos <- data[,c(1:21)]
      return(argos)
    }
  }

  ###############################################################################
  ####################run filters and perform BOD if desired#####################
  ###############################################################################

  #run MRD
  if (method == 'MRD') {
    argos <- MRD(argos, maxredun, keep_lc, keeplast, skiploc)
  } else {
    #run DAR
    if (method == 'DAR') {
      argos <- DAR(argos, maxredun, ratecoef, minrate, keep_lc)
    } else {
      #run HYB
      argos <- HYB(argos, keep_lc, maxredun,
                   keeplast, skiploc, #inputs for MRD filter
                   minrate, ratecoef, r_only, #inputs for DAR filter
                   xmigrate, xoverrun, xdirect, xangle, xpercent, test_0a, test_bz) #inputs for HYB filter
    }
  }
  outliers <- rbind(outliers, argos[which(argos$outlier == TRUE),])
  argos <- argos[which(argos$outlier == FALSE),]

  #perform best of day post-filtering if desired
  retained <- data.frame()
  if (best_of_day == TRUE) {
    stop('best-of-day filtering is only set up for data for LS filtered argos inputs, and thus not working now')
    if (is.null(minoffh)) {
      for (m in unique(lubridate::month(argos$Date))) {
        margos <- argos[which(lubridate::month(argos$Date) == m),] #all rows from the given month
        for (d in unique(lubridate::day(margos$Date))) {
          mdargos <- margos[which(lubridate::month(margos$Date) == d),] #all rows from a given day in that month
          if (nrow(mdargos) == 1) { #if only one row exists, that is the best
            retained <- rbind(retained, mdargos)
            next
          }
          if ((rankmeth != 3) & (rankmeth != 2) & (rankmeth != 1)) {
            stop('input for rankmeth must be either 1, 2, or 3')
          } else {
            if (rankmeth == 1) { #lc, iqx, iqy, nbmes
              #lc
              new_mdargos <- mdargos[which(mdargos$Quality == max(mdargos$Quality)),]
              if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                retained <- rbind(retained, new_mdargos)
                next
              }

              #iqx
              new_mdargos <- new_mdargos[which(new_mdargos$iqx == max(new_mdargos$iqx)),]
              if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                retained <- rbind(retained, new_mdargos)
                next
              }

              #iqy
              new_mdargos <- new_mdargos[which(new_mdargos$iqy == max(new_mdargos$iqy)),]
              if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                retained <- rbind(retained, new_mdargos)
                next
              }

              #nbmes
              new_mdargos <- new_mdargos[which(new_mdargos$nbmes == max(new_mdargos$nbmes)),]
              if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                retained <- rbind(retained, new_mdargos)
              }
            } else {
              if (rankmeth == 2) { #lc, iqx, nbmes, iqy
                #lc
                new_mdargos <- mdargos[which(mdargos$Quality == max(mdargos$Quality)),]
                if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                  retained <- rbind(retained, new_mdargos)
                  next
                }

                #iqx
                new_mdargos <- new_mdargos[which(new_mdargos$iqx == max(new_mdargos$iqx)),]
                if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                  retained <- rbind(retained, new_mdargos)
                  next
                }

                #nbmes
                new_mdargos <- new_mdargos[which(new_mdargos$nbmes == max(new_mdargos$nbmes)),]
                if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                  retained <- rbind(retained, new_mdargos)
                  next
                }

                #iqy
                new_mdargos <- new_mdargos[which(new_mdargos$iqy == max(new_mdargos$iqy)),]
                if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                  retained <- rbind(retained, new_mdargos)
                }
              } else {
                if (rankmeth == 3) { #lc, nbmes
                  #lc
                  new_mdargos <- mdargos[which(mdargos$Quality == max(mdargos$Quality)),]
                  if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                    retained <- rbind(retained, new_mdargos)
                    next
                  }

                  #nbmes
                  new_mdargos <- new_mdargos[which(new_mdargos$nbmes == max(new_mdargos$nbmes)),]
                  if (nrow(new_mdargos) == 1) { #if only one row exists, that is the best
                    retained <- rbind(retained, new_mdargos)
                  }
                }
              }
            }
          }
        }
      }
    } else { #use PTT duty-cycle method
      #group duty cycles into elements of a list
      dutycycles <- replicate(nrow(argos), data.frame())
      list_element <- 1
      #add first argos location to the first duty cycle list element
      dutycyles[[list_element]] <- data.frame(argos[1,])
      #loop through argos locations and place locations in respective list elements
      for (l in 1:(nrow(argos) - 1)) {
        #calculate change in time between consecutive locations
        dt <- difftime(argos$Date[(l+1)], argos$Date[l], units = 'hours')

        #if time between next location and current location is greater than the minoffh,
        #   then add the next location to the next duty cycle list element, otherwise
        #   add the next location to the current duty cycle list element
        if (dt <= minoffh) {
          dutycyles[[list_element]] <- rbind(dutycyles[[list_element]], argos[(l+1),])
        } else {
          list_element <- list_element + 1
          dutycyles[[list_element]] <- dara.frame(argos[(l+1),])
        }
      }

      #go through duty cycle list elements and find best location in each bin
      dutycyles <- dutycyles[[which(!is.null(dutycyles))]] #remove unused list elements
      retained <- data.frame()
      for (d in 1:list_element) {
        new_argos <- dutycycles[[d]]
        if (nrow(new_argos) == 1) { #if only one row exists, that is the best
          retained <- rbind(retained, new_argos)
          next
        }
        if ((rankmeth != 3) & (rankmeth != 2) & (rankmeth != 1)) {
          stop('input for rankmeth must be either 1, 2, or 3')
        } else {
          if (rankmeth == 1) { #lc, iqx, iqy, nbmes
            #lc
            new_argos <- new_argos[which(new_argos$Quality == max(new_argos$Quality)),]
            if (nrow(new_argos) == 1) { #if only one row exists, that is the best
              retained <- rbind(retained, new_argos)
              next
            }

            #iqx
            new_argos <- new_argos[which(new_argos$iqx == max(new_argos$iqx)),]
            if (nrow(new_argos) == 1) { #if only one row exists, that is the best
              retained <- rbind(retained, new_argos)
              next
            }

            #iqy
            new_argos <- new_argos[which(new_argos$iqy == max(new_argos$iqy)),]
            if (nrow(new_argos) == 1) { #if only one row exists, that is the best
              retained <- rbind(retained, new_argos)
              next
            }

            #nbmes
            new_argos <- new_argos[which(new_argos$nbmes == max(new_argos$nbmes)),]
            if (nrow(new_argos) == 1) { #if only one row exists, that is the best
              retained <- rbind(retained, new_argos)
            }
          } else {
            if (rankmeth == 2) { #lc, iqx, nbmes, iqy
              #lc
              new_argos <- new_argos[which(new_argos$Quality == max(new_argos$Quality)),]
              if (nrow(new_argos) == 1) { #if only one row exists, that is the best
                retained <- rbind(retained, new_argos)
                next
              }

              #iqx
              new_argos <- new_argos[which(new_argos$iqx == max(new_argos$iqx)),]
              if (nrow(new_argos) == 1) { #if only one row exists, that is the best
                retained <- rbind(retained, new_argos)
                next
              }

              #nbmes
              new_argos <- new_argos[which(new_argos$nbmes == max(new_argos$nbmes)),]
              if (nrow(new_argos) == 1) { #if only one row exists, that is the best
                retained <- rbind(retained, new_argos)
                next
              }

              #iqy
              new_argos <- new_argos[which(new_argos$iqy == max(new_argos$iqy)),]
              if (nrow(new_argos) == 1) { #if only one row exists, that is the best
                retained <- rbind(retained, new_argos)
              }
            } else {
              if (rankmeth == 3) { #lc, nbmes
                #lc
                new_argos <- new_argos[which(new_argos$Quality == max(new_argos$Quality)),]
                if (nrow(new_argos) == 1) { #if only one row exists, that is the best
                  retained <- rbind(retained, new_argos)
                  next
                }

                #nbmes
                new_argos <- new_argos[which(new_argos$nbmes == max(new_argos$nbmes)),]
                if (nrow(new_argos) == 1) { #if only one row exists, that is the best
                  retained <- rbind(retained, new_argos)
                }
              }
            }
          }
        }
      }
    }
  } else { #if not using one of the best_of_day methods, retained locations are all of the argos locations kept thus far
    retained <- argos
  }

  #remove outliers from retained data and create combined dataframe and make sure user locations are kept in retained
  for (i in 1:nrow(argos)) {
    if (argos$Date[i] %in% retained$Date == FALSE) {
      argos$outlier[i] <- TRUE
    }
  }
  outliers <- rbind(outliers, argos[which(argos$outlier == TRUE),])
  if (nrow(user) > 0) {
    retained <- rbind(retained, user)
  }

  #revert lc to characters
  if (nrow(retained) > 0) {
    for (r in 1:nrow(retained)) {
      if (retained$Quality[r] == 4) {
        retained$Quality[r] <- 'User'
      }
      if (retained$Quality[r] == -3) {
        retained$Quality[r] <- 'Z'
      }
      if (retained$Quality[r] == -2) {
        retained$Quality[r] <- 'B'
      }
      if (retained$Quality[r] == -1) {
        retained$Quality[r] <- 'A'
      }
    }
  }
  if (nrow(outliers) > 0) {
    for (r in 1:nrow(outliers)) {
      if (outliers$Quality[r] == -3) {
        outliers$Quality[r] <- 'Z'
      }
      if (outliers$Quality[r] == -2) {
        outliers$Quality[r] <- 'B'
      }
      if (outliers$Quality[r] == -1) {
        outliers$Quality[r] <- 'A'
      }
    }
  }
  outliers <- outliers[which(outliers$Quality != 4)]

  #create dataframe with all argos locations
  all <- rbind(outliers, retained)

  #order data by datetime column
  outliers <- outliers[order(outliers$Date),]
  retained <- retained[order(retained$Date),]
  all <- all[order(all$Date),]

  #return list of all, retained, and outliers
  return(list(all = all, retained = retained, outliers = outliers))
}
