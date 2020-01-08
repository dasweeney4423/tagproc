#' Format time data
#'
#' This function formats time data from many formats into a single standard format that is used by other tools in this package
#' @param dt The vector of times to be reformatted
#' @return A list of 10 elements:
#' \itemize{
#' \item{\strong{raw: Datetime vector NOT strptime'd.}}
#' \item{\strong{yearfirst: Datetime vector that has NOT been strptime'd, but has the year come first instead of month.}}
#' \item{\strong{strp: Datetime vector that HAS been strptime'd.}}
#' \item{\strong{mm: Vector of months.}}
#' \item{\strong{dd: Vector of days.}}
#' \item{\strong{yyyy: Vector of years.}}
#' \item{\strong{hh: Vector of hours.}}
#' \item{\strong{min: Vector of minutes.}}
#' \item{\strong{ss: Vector of seconds.}}
#' \item{\strong{j: Vector of julian days (001-366).}}
#' }
#' @note Time turner rules: Date and time have to be in the same vector. If they occur in two different columns in your data, you must concatenate them with paste() function. Date and time must be separated by a space. Date info must be delineated by a - , / , or d. Time info must be delineated with a :. If time data does not have seconds information, ":00" will be assumed for every entry.
#' @export
#' @examples #example not yet provided, sorry :(

time.turner <- function(dt) {
  # Test to see if dt is already in correct format
  dt <- as.character(dt)

  # if dt is just a string a NAs
  if ((length(unique(dt)) == 1) & (is.na(unique(dt)[1]))) {
    dtnew <- yearfirst <- dtstrp <- mm <- dd <- yyyy <- min <- ss <- hh <- j <- rep(NA, times = length(dt))
  } else {
    error.check <- tryCatch({dtnew <- dtstrp <- yearfirst <- strptime(dt, format = "%m/%d/%Y %H:%M:%S", tz = "GMT")},
                            error = function(e) {dtnew <- yearfirst <- dtstrp <- NA},
                            finally = {error = dtnew[1]})
    if (is.na(error.check[1])) {
      ##################################################
      # Assume you have unfixed width time/date data
      # Split entries into a date and a time
      dtsplit <- strsplit(dt, " ")
      split1 <- split2 <- vector()
      for (i in 1:length(dtsplit)) {
        spliti <- dtsplit[[i]]
        split1[i] <- spliti[1]
        split2[i] <- spliti[2]
      }

      # Determine which container is date, which is time
      if (length(grep(":", split1)) > 0) {
        vtime <- split1
        vdate <- split2
      } else {
        vtime <- split2
        vdate <- split1
      }

      ##################################################
      # Try to deal with time NAs
      if (length(strsplit(vtime[!is.na(vtime)], ":")[[1]]) == 2) {vtime[which(is.na(vtime))] <- "NA:NA"}
      if (length(strsplit(vtime[!is.na(vtime)], ":")[[1]]) == 3) {vtime[which(is.na(vtime))] <- "NA:NA:NA"}

      ##################################################
      # Format date
      date.goods <- which(!is.na(vdate))
      split.char <- NA
      testsplit <- strsplit(vdate[date.goods[1]], "-")[[1]]
      if (length(testsplit) > 1) {
        split.char <- "-"
      } else {
        testsplit <- strsplit(vdate[date.goods[1]], "d")[[1]]
        if (length(testsplit) > 1) {
          split.char <- "d"
        } else {
          split.char <- "/"
        }
      }
      datesplit <- strsplit(vdate, split.char)
      split1 <- split2 <- split3 <- vector()
      for (i in 1:length(datesplit)) {
        spliti <- datesplit[[i]]
        split1[i] <- spliti[1]
        split2[i] <- spliti[2]
        split3[i] <- spliti[3]
      }

      months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
      if (any(as.numeric(split1[!is.na(split1)]) > 32)) {
        yyyy <- split1
        mm <- split2
        dd <- split3
      } else {
        yyyy <- split3
        if (split2[date.goods[1]] %in% months) {
          mm <- split2
          dd <- split1
        } else {
          mm <- split1
          dd <- split2
        }
      }
      mm[mm==months[1]] <- "01"
      mm[mm==months[2]] <- "02"
      mm[mm==months[3]] <- "03"
      mm[mm==months[4]] <- "04"
      mm[mm==months[5]] <- "05"
      mm[mm==months[6]] <- "06"
      mm[mm==months[7]] <- "07"
      mm[mm==months[8]] <- "08"
      mm[mm==months[9]] <- "09"
      mm[mm==months[10]] <- "10"
      mm[mm==months[11]] <- "11"
      mm[mm==months[12]] <- "12"
      mm[nchar(mm) < 2 & !is.na(mm)] <- paste0("0", mm[nchar(mm) < 2 & !is.na(mm)])
      dd[nchar(dd) < 2 & !is.na(dd)] <- paste0("0", dd[nchar(dd) < 2 & !is.na(dd)])

      for (Y in 1:length(yyyy)) {
        if ((nchar(yyyy[Y]) < 3) & (!is.na(yyyy[Y]))) {
          if (as.numeric(yyyy[Y]) > 50) {
            yyyy[Y] <- paste0("19", yyyy[Y])
          } else {
            yyyy[Y] <- paste0("20", yyyy[Y])
          }
        }
      }

      ##################################################
      # Format time
      timesplit <- strsplit(vtime, ":")
      split1 <- split2 <- split3 <- vector()
      for (i in 1:length(timesplit)) {
        spliti <- timesplit[[i]]
        split1[i] <- spliti[1]
        split2[i] <- spliti[2]
        if (length(spliti) > 2) {split3[i] <- spliti[3]}
      }

      hh  <- split1
      hh[nchar(hh) < 2] <- paste0("0", hh[nchar(hh) < 2])

      min <- split2
      min[nchar(min) < 2] <- paste0("0", min[nchar(min) < 2])

      ss <- split3
      if (length(ss) == 0) {
        ss <- rep("00", times = length(min))
      } else {
        ss[nchar(ss) < 2] <- paste0("0", ss[nchar(ss) < 2])
      }

      ##################################################
      # Concatenate date and time, then change class
      dtnew <- paste0(mm, "/", dd, "/", yyyy, " ", hh, ":", min, ":", ss)
      yearfirst <- paste0(yyyy, "/", mm, "/", dd, " ", hh, ":", min, ":", ss)
      dtstrp <- strptime(dtnew, format = "%m/%d/%Y %H:%M:%S", tz = "GMT")

      # Calculate julian day
      j <- strftime(dtstrp, format = "%j")
    } else { # end of IF error=="error" check
      yyyy <- substr(dtnew, 1, 4)
      mm <- substr(dtnew, 6, 7)
      dd <- substr(dtnew, 9, 10)
      hh <- substr(dtnew, 12, 13)
      min <- substr(dtnew, 15, 16)
      ss <- substr(dtnew, 18, 19)
      j <- strftime(dtstrp, format = "%j")
    } # end of else of IF error=="error" check
  } # end of else for if dt is all NA

  # Return a list: first element is the vector NOT strptime'd. The second element HAS been strptime'd
  return(list(raw = dtnew,
              yearfirst = yearfirst,
              strp = dtstrp,
              mm = mm,
              dd = dd,
              yyyy = yyyy,
              hh = hh,
              min = min,
              ss = ss,
              j = j))
}
