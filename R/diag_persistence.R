#'Test for functional group persistence
#'
#' Inspects each time point. If at any point in time a species falls below a predefined floor it
#' is flagged. The floors are specifies as a proportion of initial biomass. The last n years of the run are
#' used in the test
#'
#'@param fgs A character string. Path to location of functional groups file.
#'@param biomind A character string. Path to the BiomIndx.txt file.
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for persistence.
#'(Default = NULL, uses all species with \code{IsTurnedOn=1} in \code{fgs} file)
#'@param nYrs Numeric scalar. Number of years from the end of the time series that persistence must occur.
#' (Default = NULL, persistence must occur throughout entire time series)
#'@param floor Numeric scalar. Proportion of initial biomass for which for which all species
#'are measured against. (Default = 0.1, all species need to be above 10% of initial biomass). Range should be 0-1
#'@param display Boolean. Flag to indicate whether to return only species that pass the test, fail the test, or all species (Default = NULL)
#'@param tol Numeric scalar. Tolerance level to add to biomass floor (Default = 1E-6)
#'
#'
#'@return Returns a data frame of species which do not meet defined persistence criteria.
#'
#'\item{species}{The common name of the species/functional group}
#'\item{code}{Atlantis Code for species/functional group}
#'\item{initialBiomass}{Starting value of Biomass for species/functional group}
#'\item{minimumBiomass}{The smallest value of biomass observed in the run}
#'\item{tminimumBiomass}{The time step in which \code{minimumBiomass} occurred }
#'\item{proportionInitBio }{The proportion of initial biomass at \code{tminimumBiomass}}
#'\item{t1}{The first time step persistence was not met}
#'\item{tn}{The last time step persistence was not met}
#'\item{nts}{The total number of time steps that persistence was not met}
#'\item{pass}{Boolean. Indicates if the Atlantis group passed persistence test}
#'
#'@export
#'
#'@family diags
#'
#'@importFrom magrittr %>%
#'@importFrom rlang .data
#'
#'@examples
#'\dontrun{
#'
#'# Declare paths to files required
#' biomind <- paste("Full path to file","xxxBiomIndx.txt")
#' fgs <- paste("Full path to file","functioalGgroups.csv")
#'
#' # find all species that do not have biomass > 0 for any time during the run.
#' diag_persistence(fgs,biomind,speciesCodes=NULL, nYrs = NULL, floor = 0)
#'
#' # only evaluate herring. Require stability over the last 10 years of the run and all values should
#' # exceed 10% of initial biomass
#' diag_persistence(fgs,biomind, speciesCodes="HER", nYrs = 10, floor = 0.1)
#'
#'}

diag_persistence <- function(fgs,
                             biomind,
                             speciesCodes=NULL,
                             nYrs = NULL,
                             floor = 0.1,
                             display=NULL,
                             tol = 1E-6){


  # read in biomass data and qualify species Codes
  biom <- get_model_biomass(fgs,biomind,speciesCodes)
  modelBiomass <- biom$modelBiomass
  speciesCodes <- biom$speciesCodes

  # find final time step value
  maxRuntime <- max(modelBiomass$time)

  if (is.null(nYrs)) { # use all time series
    filterTime <- 0
  } else { # last n years
    filterTime <- maxRuntime - (365*nYrs)
  }

  # For each species calculate which time steps biomass is below persistence threshold
  # we look at biomass < % initial Biomass
  status <- modelBiomass %>%
    dplyr::filter(.data$code %in% speciesCodes) %>%
    dplyr::select(.data$code,.data$species, .data$time, .data$atoutput) %>%
    dplyr::group_by(.data$code) %>%
    dplyr::mutate(initialBiomass = dplyr::first(.data$atoutput)) %>%
    dplyr::filter(.data$time >= filterTime) %>%
    dplyr::mutate(proportionInitBio = dplyr::if_else(is.nan(.data$atoutput/.data$initialBiomass),0,.data$atoutput/.data$initialBiomass)) %>%
    dplyr::mutate(proportionInitBio = as.numeric(trimws(format(round(.data$proportionInitBio,3),nsmall=3)))) %>%
    #    dplyr::filter(proportionInitBio <= (floor + tol)) %>%
    dplyr::mutate(pass = .data$proportionInitBio >= (floor + tol)) %>%
    dplyr::ungroup()

  # num times threshold exceeded, when largest exceedance occurs and value of biomass,
  # range of exceedances

  persistenceF <- status %>%
    dplyr::filter(.data$pass == F) %>%
    dplyr::group_by(.data$code) %>%
    dplyr::mutate(minimumBiomass = min(.data$atoutput)) %>%
    dplyr::mutate(nts = dplyr::n()) %>%
    dplyr::mutate(tminimumBiomass = dplyr::case_when(.data$atoutput == min(.data$atoutput) ~ .data$time)) %>%
    dplyr::mutate(t1 = min(.data$time), tn = max(.data$time)) %>%
    dplyr::filter(!is.na(.data$tminimumBiomass)) %>%
    dplyr::select(.data$code,.data$species,.data$initialBiomass,
                  .data$proportionInitBio, .data$minimumBiomass,
                  .data$tminimumBiomass, .data$t1, .data$tn, .data$nts,.data$pass) %>%
    dplyr::filter(.data$tminimumBiomass == min(.data$tminimumBiomass)) %>%
    dplyr::ungroup()

  codesFailed <- persistenceF %>%
    dplyr::pull(.data$code)

  persistenceT <- status %>%
    dplyr::filter(!(.data$code %in% codesFailed)) %>%
    dplyr::group_by(.data$code) %>%
    dplyr::mutate(minimumBiomass = min(.data$atoutput)) %>%
    dplyr::mutate(nts = 0) %>%
    dplyr::mutate(tminimumBiomass = dplyr::case_when(.data$atoutput == min(.data$atoutput) ~ .data$time)) %>%
    dplyr::mutate(t1 = NA, tn = NA) %>%
    dplyr::filter(!is.na(.data$tminimumBiomass)) %>%
    dplyr::select(.data$code,.data$species,.data$initialBiomass,.data$proportionInitBio,
                  .data$minimumBiomass, .data$tminimumBiomass, .data$t1, .data$tn, .data$nts,.data$pass) %>%
    dplyr::filter(.data$tminimumBiomass == min(.data$tminimumBiomass)) %>%
    dplyr::ungroup()

  persistence <- rbind(persistenceT,persistenceF)


  if (is.null(display)) {
    # return all
  } else if (display) {
    # return species that pass
    persistence <- persistence %>% dplyr::filter(.data$pass==T)

  } else {
    # return all that fail
    persistence <- persistence %>% dplyr::filter(.data$pass!=T)

  }

  persistence <- persistence %>%
    dplyr::arrange(.data$pass,.data$proportionInitBio)

  return(persistence)

}
