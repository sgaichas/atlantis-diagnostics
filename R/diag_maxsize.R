#' Test for max fish size
#'
#' Calculate average size of fish through time and check to see if any fish exceed the
#' reported largest size.
#'
#'@param atlDir A character string. Path to location of model run output files.
#'@param runPrefix A character string. Prefix added to all atlantis output files.
#'@param speciesStats Data frame.
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for persistence.
#'(Default = NULL, uses all species)
#'@param nYrs Numeric scalar. Number of years from the end of the time series that persistence must occur.
#' (Default = NULL, persistence must occur throughout entire time series)
#'
#'@importFrom magrittr %>%
#'@importFrom rlang .data
#'
#'
#'@return Returns a data frame of species which do not meet defined max size criterion.
#'
#'\item{species}{The common name of the species/functional group}
#'\item{code}{Atlantis Code for species/functional group}

#'@export
#'
#'
#'@examples
#'\dontrun{
#'
#'}

diag_maxsize <- function(atlDir,runPrefix,speciesStats,speciesCodes=NULL, nYrs = NULL){

  # look for AnnualAgeBiomassIndx.txt and AnnualAgeNumbersInd.txt in output folder
  biomassFile <- paste0(atlDir,runPrefix,"AnnualAgeBiomIndx.txt")
  ageFile <- paste0(atlDir,runPrefix,"AnnualAgeNumbersIndx.txt")
  if (!file.exists(biomassFile)) {
    stop(paste("Can't find Biomass file - ",biomassFile))
  }
  if (!file.exists(ageFile)) {
    stop(paste("Can't find age file - ",ageFile))
  }

  biomData <- atlantistools::load_txt(biomassFile) %>%
    atlantistools::preprocess_txt(.,into = "code",removeZeros=F) %>%
    dplyr::filter(code %in% speciesCodes) %>%
    dplyr::rename(biomass = .data$atoutput)

  ageData <- atlantistools::load_txt(ageFile) %>%
    atlantistools::preprocess_txt(.,into = "code",removeZeros=F) %>%
    dplyr::filter(.data$code %in% speciesCodes) %>%
    dplyr::rename(age = .data$atoutput)


  joinData <- dplyr::left_join(biomData,ageData,by = c("time","code")) %>%
    dplyr::mutate(meanWeight=.data$biomass/.data$age) %>%
    dplyr::group_by(code) %>%
    dplyr::summarise(maxMeanWeight = max(meanWeight))

  # filter by speciesCodes

  # then by data availability

  # compare to stats

  # report diagnostic

  return()

}




  # find final time step value
  maxRuntime <- max(modelBiomass$time)

  # check for valid species Codes & clean
  speciesCodes <- check_species_codes(modelBiomass,speciesCodes)

  if (is.null(nYrs)) { # use all time series
    filterTime <- 0
  } else { # last n years
    filterTime <- maxRuntime - (365*nYrs)
  }

  # For each species calculate which time steps biomass is below persistence threshold
  # we look at biomass < % initial Biomass
  status <- modelBiomass %>%
    dplyr::filter(code %in% speciesCodes) %>%
    dplyr::select(code,species, time, atoutput) %>%
    dplyr::group_by(code) %>%
    dplyr::mutate(initialBiomass = dplyr::first(atoutput)) %>%
    dplyr::filter(time >= filterTime) %>%
    dplyr::mutate(proportionInitBio = dplyr::if_else(is.nan(atoutput/initialBiomass),0,atoutput/initialBiomass)) %>%
    dplyr::mutate(proportionInitBio = as.numeric(trimws(format(round(proportionInitBio,3),nsmall=3)))) %>%
    #    dplyr::filter(proportionInitBio <= (floor + tol)) %>%
    dplyr::mutate(pass = proportionInitBio >= (floor + tol)) %>%
    dplyr::ungroup()

  # num times threshold exceeded, when largest exceedance occurs and value of biomass,
  # range of exceedances

  persistenceF <- status %>%
    dplyr::filter(pass == F) %>%
    dplyr::group_by(code) %>%
    dplyr::mutate(minimumBiomass = min(atoutput)) %>%
    dplyr::mutate(nts = dplyr::n()) %>%
    dplyr::mutate(tminimumBiomass = dplyr::case_when(atoutput == min(atoutput) ~ time)) %>%
    dplyr::mutate(t1 = min(time), tn = max(time)) %>%
    dplyr::filter(!is.na(tminimumBiomass)) %>%
    dplyr::select(code,species,initialBiomass,proportionInitBio, minimumBiomass, tminimumBiomass, t1, tn, nts,pass) %>%
    dplyr::filter(tminimumBiomass == min(tminimumBiomass)) %>%
    dplyr::ungroup()

  codesFailed <- persistenceF %>%
    dplyr::pull(code)

  persistenceT <- status %>%
    dplyr::filter(!(code %in% codesFailed)) %>%
    dplyr::group_by(code) %>%
    dplyr::mutate(minimumBiomass = min(atoutput)) %>%
    dplyr::mutate(nts = 0) %>%
    dplyr::mutate(tminimumBiomass = dplyr::case_when(atoutput == min(atoutput) ~ time)) %>%
    dplyr::mutate(t1 = NA, tn = NA) %>%
    dplyr::filter(!is.na(tminimumBiomass)) %>%
    dplyr::select(code,species,initialBiomass,proportionInitBio, minimumBiomass, tminimumBiomass, t1, tn, nts,pass) %>%
    dplyr::filter(tminimumBiomass == min(tminimumBiomass)) %>%
    dplyr::ungroup()

  persistence <- rbind(persistenceT,persistenceF)


  if (is.null(display)) {
    # return all
  } else if (display) {
    # return species that pass
    persistence <- persistence %>% dplyr::filter(pass==T)

  } else {
    # return all that fail
    persistence <- persistence %>% dplyr::filter(pass!=T)

  }

  persistence <- persistence %>%
    dplyr::arrange(pass,proportionInitBio)

  return(persistence)

}
