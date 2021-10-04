#'Test for functional group stability
#'
#'\code{diag_stability} determines whether the model reaches a steady state
#'over the last n years of a run. Stability is loosely defined as a species/groups biomass reaching a
#'stable level (measured as having tolerable trend in a defined time range).
#'
#'
#'@param fgs A character string. Path to location of functional groups file.
#'@param biomind A character string. Path to the BiomIndx.txt file.
#'@param initialYr Numeric Scalar. Year in which the model run was initiated. (Default = 1964)
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for stability.
#'(Default = NULL, uses all species found in  \code{modelBiomass})
#'@param nYrs Numeric scalar. Number of years from the end of the time series that stability must occur.
#' (Default = 20 years)
#'@param relChangeThreshold Numeric Scalar. Maximum magnitude of relative change of slope (Default = 0.01)
#'
#'@return Returns a data frame of all species and how they measure up against the stability criterion
#'\item{code}{Atlantis Code for species/functional group}
#'\item{species}{The common name of the species/functional group}
#'\item{t1Fit}{Value of fitted biomass for the first of year data used in the fit}
#'\item{mtPerYear}{Double. The value of the slope parameter (year) }
#'\item{relChange}{Rate of increase relative to \code{t1Fit}(\code{mtPerYear}/\code{t1Fit})}
#'\item{aveBio}{mean biomass for the last \code{nYrs} years }
#'\item{pass}{Logical. Does the species/group pass the test for stability}
#'
#'
#'@section Details:
#'
#'Formally the following model is fit to the last n years of the run:
#'
#' \deqn{biomass_t = \mu + \beta.t + \epsilon_t  where \epsilon_t ~ IID N(0,\sigma^2)}
#'
#' where null hypothesis,  \deqn{H0:\beta=0}
#'
#' Note: annual biomass is used in fitting. Species with mean annual biomass < 1 metric ton over the last
#' n years of the run are not considered stable. They are reported to Fail the test and NaNs returned
#'
#'
#'@family diags
#'
#'@export
#'
#'@importFrom magrittr %>%
#'@importFrom rlang .data
#'
#'@examples
#'\dontrun{
#'# Declare paths to files required
#'
#' biomind <- paste("Full path to file","xxxBiomIndx.txt")
#' fgs <- paste("Full path to file","functioalGgroups.csv")
#'
#' # Perform stability test on all species/groups using the last 20 years of the run
#' diag_stability(fgs, biomind, nYrs = 20)
#'
#' # Only perform test on herring and white hake.
#' # Require stability over the last 10 years of the run and and use a
#' # relative change in the slope = 0.01 as the criterion for pass or fail
#' diag_stability(fgs,biomind, speciesCodes=c("HER","WHK"), nYrs = 10, relChangeThreshold = 0.01)
#'}

diag_stability <- function(fgs,
                           biomind,
                           initialYr = 1964,
                           speciesCodes,
                           nYrs = 20,
                           relChangeThreshold = 0.01){

  # read in biomass data and qualify species Codes
  biom <- get_model_biomass(fgs,biomind,speciesCodes)
  modelBiomass <- biom$modelBiomass
  speciesCodes <- biom$speciesCodes

  # determine time frame in which to perform "test"
  # filter model biomass and average over the year
  # join with initial biomass
  modelBiomass <- modelBiomass %>%
    dplyr::select(.data$time, .data$atoutput, .data$code, .data$species) %>%
    dplyr::rename(value=.data$atoutput) %>%
    dplyr::mutate(yearStep= floor(.data$time/365)) %>%
    dplyr::group_by(.data$code,.data$yearStep,.data$species) %>%
    dplyr::summarise(meanBio = mean(.data$value),.groups="drop") %>% # average over year, removes seasonal cycle
    dplyr::mutate(year = .data$yearStep+initialYr, yearStep=NULL) %>%
    dplyr::ungroup()

  if (is.null(nYrs)) { # use all time series
    filterTime <- initialYr
  } else { # last n years
    filterTime <- max(modelBiomass$year) - nYrs
  }

  # select species, time frame, group, nest
  stable <- modelBiomass %>%
    dplyr::filter(.data$year > filterTime) %>%
    dplyr::filter(.data$code %in% speciesCodes) %>%
    dplyr::group_by(.data$code,.data$species) %>%
    #dplyr::mutate(t1Biomass = dplyr::first(meanBio)) %>%
    dplyr::group_by(.data$code,.data$species) %>%
    tidyr::nest() # produces a column called data

  # fit linear model (biomass = alpha + year.beta) to each species over the time frame specified
  # then extract slope and significance from model fit

  stability <- stable %>%
    dplyr::mutate(model = purrr::map(.data$data,fitlm)) %>%
    #dplyr::transmute(species,mtPerYear = purrr::map_dbl(model,coefs),pValue=purrr::map_dbl(model,pVals)) %>%
    dplyr::transmute(.data$species,t1Fit = purrr::map_dbl(.data$model,fittedVal),
                     mtPerYear = purrr::map_dbl(.data$model,coefs),
                     aveBio = purrr::map_dbl(.data$data,meanData)) %>%
    #dplyr::mutate(pass = pValue > sigTest) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(relChange = .data$mtPerYear/.data$t1Fit) %>%
    dplyr::mutate(sign = base::sign(.data$mtPerYear)) %>%
    dplyr::mutate(pass = abs(.data$relChange) < relChangeThreshold) %>%
    # average biomass < 1 metric ton
    dplyr::mutate(pass = dplyr::if_else(.data$aveBio < 1,F,.data$pass)) %>%
    dplyr::mutate(sign = dplyr::if_else(.data$aveBio < 1,NaN,.data$sign)) %>%
    dplyr::mutate(relChange = dplyr::if_else(.data$aveBio < 1,NaN,.data$relChange)) %>%
    dplyr::mutate(t1Fit = dplyr::if_else(.data$aveBio < 1,NaN,.data$t1Fit)) %>%
    dplyr::mutate(mtPerYear = dplyr::if_else(.data$aveBio < 1,NaN,.data$mtPerYear))


  stability <- stability %>%
    dplyr::arrange(.data$pass,dplyr::desc(abs(.data$relChange)))

  return(stability)



}









