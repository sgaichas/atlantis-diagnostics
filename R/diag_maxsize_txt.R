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
#'
#'@noRd
#'
#'
#'@examples
#'\dontrun{
#'
#'}

diag_maxsize_txt <- function(atlDir,runPrefix,speciesStats,speciesCodes=NULL, nYrs = NULL){

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
    atlantistools::preprocess_txt(.data,into = "code",removeZeros=F) %>%
    dplyr::filter(.data$code %in% speciesCodes) %>%
    dplyr::rename(biomass = .data$atoutput)

  ageData <- atlantistools::load_txt(ageFile) %>%
    atlantistools::preprocess_txt(.data,into = "code",removeZeros=F) %>%
    dplyr::filter(.data$code %in% speciesCodes) %>%
    dplyr::rename(age = .data$atoutput)


  joinData <- dplyr::left_join(biomData,ageData,by = c("time","code")) %>%
    dplyr::mutate(meanWeight=.data$biomass/.data$age) %>%
    dplyr::group_by(.data$code) %>%
    dplyr::summarise(maxMeanWeight = max(.data$meanWeight))

  # filter by speciesCodes

  # then by data availability

  # compare to stats

  # report diagnostic for species have max
  # report largest size of all
  # report no data for others

  return()

}

