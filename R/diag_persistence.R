#'Test for functional group persistence
#'
#' Inspects each time point. If at any point in time a species falls below a predefined floor it
#' is flagged. The floors are specifies as a proportion of initial biomass. The last n years of the run are
#' used in the test
#'
#'@param modelBiomass A data frame. Total biomass of all groups over time, read in from
#'Atlantis ...BioInd.txt output using \code{atlantisom::load_bioind}.
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for persistence.
#'(Default = NULL, uses all species found in  \code{modelBiomass})
#'@param nYrs Numeric scalar. Number of years from the end of the time series that persistence must occur.
#' (Default = NULL, persistence must occur throughout entire time series)
#'@param floor Numeric scalar. Proportion of initial biomass for which for which all species
#'are measured against. (Default = 0, all species need to be non zero). Range should be 0-1
#'@param tol Numeric scalar. Tolerance level to add to biomass floor (Default = 1E-6)
#'@param plot A logical value specifying if the function should generate plots or
#'not. (Default = F).
#'
#'@importFrom magrittr %>%
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
#'
#'@export
#'
#'
#'@examples
#'\dontrun{
#'# Declare paths to files required
#' biol.file <- "neus_outputBiomIndx.txt"
#' file_fgs <- "neus_groups.csv"
#' # use atlantisom to read them in
#' fgs <- atlantisom::load_fgs(inDir,file_fgs)
#' modelBiomass <- atlantisom::load_bioind(outDir,biol.file,fgs)
#'
#' # find all species that do not have biomass > 0 for any time during the run.
#'
#' diag_persistence(modelBiomass,speciesCodes=NULL, nYrs = NULL, floor = 0)
#'
#' # only evaluate herring. Require stability over the last 10 years of the run and all values should
#' exceed 10% of initial biomass
#' diag_persistence(modlBiomass, speciesCodes="HER", nYrs = 10, floor = 0.1)
#'
#'}

diag_persistence <- function(modelBiomass, speciesCodes=NULL, nYrs = NULL, floor = 0, tol = 1E-6, plot=F){

  # need in annual units? Or fail when any output timestep below threshold?
  # make safe for migratory species, assume that over the course of the year mean B > 0.
  # assumes biomass never goes negative in atlantis

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
    dplyr::filter(time >= filterTime) %>%
    dplyr::filter(code %in% speciesCodes) %>%
    dplyr::select(code,species, time, atoutput) %>%
    dplyr::group_by(code) %>%
    dplyr::mutate(initialBiomass = dplyr::first(atoutput)) %>%
    dplyr::mutate(proportionInitBio = atoutput/initialBiomass) %>%
    dplyr::filter(proportionInitBio <= (floor + tol)) %>%
    dplyr::ungroup()

  if (nrow(status) == 0) { # all pass species persist at desired levels
    return(persistence=NULL)
  }

  # num times threshold exceeded, when largest exceedance occurs and value of biomass,
  # range of exceedances
  persistence <- status %>%
    dplyr::group_by(code) %>%
    dplyr::mutate(minimumBiomass = min(atoutput)) %>%
    dplyr::mutate(nts = dplyr::n()) %>%
    dplyr::mutate(tminimumBiomass = dplyr::case_when(atoutput == min(atoutput) ~ time)) %>%
    dplyr::mutate(t1 = min(time), tn = max(time)) %>%
    dplyr::filter(!is.na(tminimumBiomass)) %>%
    dplyr::select(code,species,initialBiomass,proportionInitBio, minimumBiomass, tminimumBiomass, t1, tn, nts) %>%
    dplyr::filter(tminimumBiomass == min(tminimumBiomass)) %>%
    dplyr::ungroup()



  return(persistence)

}
