data <- read_csv("C:/Users/marec/Docs/ZcTag059-Portal-164613/WC-pretty/164613-Locations.csv")

crawl.apply <- function(data, model.interval = 1, land.adjust = TRUE, img.path = NULL) {
  if (nrow(data) > 0) {
    #convert Date to POSIXct if necessary
    if ('POSIXct' %in% class(data$Date) == FALSE) {
      data$Date <- as.POSIXct(as.character(data$Date), format = '%m/%d/%Y %H:%M:%S', tz='UTC')
    }

    ## Preparing Input Data for `crawl`
    #### Converting to spatial data
    sf_locs <- st_as_sf(data, coords = c("Longitude", "Latitude")) %>%
      st_set_crs("+init=epsg:4326")

    ### Error Parameters for GPS Data
    # set code to eliminate records with NA for error_radius. Does not eliminate any records for test tracks used.
    data <- data %>%
      dplyr::filter(!(is.na(Error.radius)))

    ### Duplicate Times
    make_unique <- function(x) {
      xts::make.time.unique(x$Date, eps = 1)
    }

    data <- data %>% dplyr::arrange(DeployID, Date)
    data$Date <- suppressWarnings(make_unique(data))

    #########################################################
    ### Course Speed Filter
    # Note that `argosfilter::sdafilter` additionally filters out very sharp turning angles.
    speed_filt <- function(x) {
      argosfilter::sdafilter(
        lat = x$Latitude,
        lon = x$Longitude,
        dtime = x$Date,
        lc = x$LocationQuality,
        vmax = 7.5,
        ang = -1
      )
    }

    data <- data %>%
      dplyr::group_by(DeployID) %>% tidyr::nest() %>%
      dplyr::mutate(filtered = purrr::map(data, speed_filt)) %>%
      tidyr::unnest(cols = c(data, filtered)) %>%
      dplyr::filter(filtered %in% c("not", "end_location")) %>%
      dplyr::select(-filtered) %>%
      dplyr::arrange(DeployID, Date)

    #########################################################
    ### Create a Spatial Object
    sf_locs <- sf::st_as_sf(data, coords = c("Longitude","Latitude")) %>%
      sf::st_set_crs(., 4326)
    # Best projection for Southern California region is abbreviated by the EPSG code of 2230.
    sf_locs <- sf::st_transform(sf_locs, 2230)
    sf::st_geometry(sf_locs)

    #At this point, `sf_locs` is a valid format for input to `crawl::crwMLE()`.
    #But run this next bit for compatibility with code further down,
    #to translate the simple feature geometry for our point coordinates into an `x` and `y` column (note the use of a local
    # function `sfc_as_cols()` to facilitate this). (This step may no longer be necessary if adjusted for later in code).

    sfc_as_cols <- function(x, names = c("x", "y")) {
      stopifnot(inherits(x, "sf") && inherits(sf::st_geometry(x), "sfc_POINT"))
      ret <- do.call(rbind, sf::st_geometry(x))
      ret <- tibble::as_tibble(ret)
      stopifnot(length(names) == ncol(ret))
      ret <- setNames(ret, names)
      dplyr::bind_cols(x, ret)
    }

    # run this line instead of snippet in crawl-practical that includes activity data
    sf_locs <- suppressWarnings(sf_locs %>% sfc_as_cols())

    #########################################################
    ## Determining Your Model Parameters
    ### Create a Nested Data Structure
    sf_locs <- sf_locs %>%
      dplyr::group_by(DeployID) %>% dplyr::arrange(Date) %>%
      tidyr::nest() %>%
      dplyr::mutate(data = purrr::map(data, sf::st_as_sf))

    ### Create model.matrix for Ellipse Errors
    ellipse_matrix <- function(x) {
      if (inherits(x, "sf")) {
        sf::st_geometry(x) <- NULL
      }
      ret <- model.matrix(
        ~ Error.Semi.major.axis + Error.Semi.minor.axis +
          Error.Ellipse.orientation,
        model.frame(~ ., x, na.action = na.pass))[,-1]
    }

    sf_locs <- sf_locs %>%
      dplyr::mutate(diag = purrr::map(data, ellipse_matrix),
                    diag = purrr::map(diag, ~ crawl::argosDiag2Cov(
                      .x[,1], .x[,2], .x[,3])),
                    data = purrr::map2(data, diag, dplyr::bind_cols)
      ) %>%
      dplyr::select(-diag)

    ### Create Model Parameters
    init_params <- function(d) {
      ret <- list(a = c(d$x[1], 0,
                        d$y[1], 0),
                  P = diag(c(10 ^ 2, 10 ^ 2,
                             10 ^ 2, 10 ^ 2)))
    }

    # Like in the book example, the first two values in the vector are set to 1 as we are providing the
    # error structure with the ellipse information. The third value in the vector is the _sigma_ and
    # this is almost always set to _NA_. The fourth value is _beta_ which is the autocorrelation parameter.
    # This is also, typically, set to _NA_ as we want to estimate this. However, in some cases, the
    # movement data will be more representative of Brownian motion and the model will be more likely
    # to fit successfully if this is fixed to a value of _4_.

    # Constraints: in this case, we let the _sigma_ parameter ranged from -Inf to Inf while
    # limiting _beta_ to range between -4 and 4.

    sf_locs <- sf_locs %>%
      dplyr::mutate(init = purrr::map(data, init_params),
                    fixpar = rep(
                      list(c(1, 1, NA, NA)),
                      nrow(.)
                    ),
                    constr = rep(list(
                      list(lower = c(-Inf, -4), upper = (c(Inf, 4)))),
                      nrow(.))
      )

    #########################################################
    ## Fitting with `crawl::crwMLE()`

    ### Create wrapper function for `crawl::crwMLE()` to align with _tidy_ approach
    # We pass our observed data in as `d` and then the parameters `init`, `fixpar`, and `constr`.
    # Additionally, we can specify `tryBrownian = FALSE` if we don't want the model
    # to fit with Brownian motion as the final try.

    # Each of the prior functions is specified within the function and the function
    # cycles through a for loop of calls to `crawl::crwMLE()` with each prior. Note,
    # the first prior value is _NULL_ to specify our first model fit try is without
    # any prior.

    # Our wrapper function also checks the observed data, `d`, for any _activity_
    # column to determine whether activity parameter should be included in themodel fit.

    fit_crawl <- function(d, init, fixpar, constr,
                          tryBrownian = TRUE) {

      priors <- list(
        NULL,
        ln_prior = function(par) { dnorm(par[2], 4, 4, log = TRUE) },
        lap_prior = function(par) { -abs(par[2] - 4) / 5 },
        reg_prior = function(par) { dt(par[2] - 3, df = 1, log = TRUE) }
      )
      #cycle through 4 different prior values. the first being no prior/NULL
      for (prior in priors) {
        fit <- crawl::crwMLE(
          mov.model =  ~ 1,
          if (any(colnames(d) == "activity")) {
            activity <- ~ I(activity)
          } else {activity <- NULL},
          err.model = list(
            x =  ~ ln.sd.x - 1,
            y =  ~ ln.sd.y - 1,
            rho =  ~ error.corr
          ),
          data = d,
          Time.name = "Date",
          initial.state = init,
          fixPar = fixpar,
          prior = prior,
          constr = constr,
          attempts = 1,
          control = list(
            trace = 0
          ),
          initialSANN = list(
            maxit = 1500,
            trace = 0
          )
        )
        if (!any(is.nan(fit$se))) {
          return(fit)
        }
      }

      #########################################################
      # at this point, the most likely reason for failing to fit is b/c the mov't is
      # more Brownian in nature. Here, we fix beta at 4 which specifies Brownian
      if (any(is.nan(fit$se)) && tryBrownian) {
        fixPar = c(1, 1, NA, 4)

        fit <- crawl::crwMLE(
          mov.model =  ~ 1,
          err.model = list(
            x =  ~ ln.sd.x - 1,
            y =  ~ ln.sd.y - 1,
            rho =  ~ error.corr
          ),
          data = d,
          Time.name = "Date",
          initial.state = init,
          fixPar = fixPar,
          attempts = 1,
          control = list(
            trace = 0
          ),
          initialSANN = list(
            maxit = 500,
            trace = 0
          )
        )
      }
    }

    #########################################################
    #Store the model results and estimated parameters as a list-column within our _nested_ tibble.
    #This section relies on the `crawl::tidy_crwFit()` function which is only available with the latest `devel` version of `crawl` from GitHub.

    data_fit <- suppressWarnings(sf_locs %>%
      dplyr::mutate(fit = purrr::pmap(list(d = data, init = init,
                                           fixpar = fixpar, constr = constr),
                                      fit_crawl),
                    params = purrr::map(fit, crawl::tidy_crwFit)))

    # Explore the resulting parameter estimates. Note that for Ziphius, proximity of beta to 4
    # means it is essentially modeled as Brownian motion.
    # This is somewhat less extreme if use angle filter in speedfilter step in preprocessing above,
    # but seems like that would discard legitimate positions in the case of these tracks (?).
    # Note that get NAs for x and y SE because fixed these at 1 since using ellipse information on uncertainty.
    # (Note: if just want to look at fit for one track, can use data_fit$fit[[1]])

    #########################################################
    ## Exploring and Troubleshooting Model Results
    ## Predicting a Movement Track
    # The `crawl::crwPredict()` function requires two key arguments: an `object.crwFit` and a `predTime` argument.
    # Note that first point of second track requires further investigation based on following diagnostic plots.

    data_fit <- data_fit %>%
      dplyr::mutate(predict = purrr::map(fit,
                                         crawl::crwPredict,
                                         predTime = paste0(model.interval, ' hour')))
    if (!is.null(img.path)) {
      pdf(file = paste0(img.path, "_", Sys.Date(), ".pdf"), height = 8, width = 7)
      data_fit$predict %>% purrr::walk(crawl::crwPredictPlot, plotType = "map")
      dev.off()
    }

    #########################################################
    ### Additional processing for spatial analysis
    # Perform additional processing of the prediction object to prepare for analysis or visualization
    # by creating a custom function `as.sf()` which will convert our `crwPredict` object
    # into an `sf` object of either "POINT" or "MULTILINESTRING"

  } # end of if nrow(data) > 0
}
