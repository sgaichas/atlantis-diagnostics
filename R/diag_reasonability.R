#'Test for all functional groups to persist within reasonable bounds
#'
#' Inspects each time point. If at any point in time a species falls below or rises above predefined
#' bounds it is flagged. The term "reasonable" is based on survey data when available.
#'
#'@param modelBiomass A data frame. Total biomass of all groups over time, read in from
#'Atlantis ...BioInd.txt output using \code{atlantisom::load_bioind}.
#'@param initialYr Numeric Scalar. Year in which the model run was initiated. (Default = 1964)
#'@param speciesNames Character vector. A vector of species names in which to test for reasonableness
#'(Default = NULL, uses all species found in  \code{biomass}. Results will depend on availability of calculable
#'bounds based on "data")
#'@param realBiomass A data frame. biomass time series (from assessments, stock SMART or otherwise) for species.
#' \code{realBiomass} should be in long format with at least 3 columns labeled Year, Value, CommonName.
#'@param nYrs Numeric scalar. Number of years from the end of the time series for which reasonableness is checked.
#' (Default = NULL, entire time series is used)
#'@param bounds Numeric vector. Size of 1x2 containing the values in which to multiple lower and upper bounds of observed data.
#'For example (Default = c(1,1)) indicating use of min and max of observed biomass
#'Proportion of initial biomass for which for which all species
#'are measured against. (Default = 0, all species need to be non zero). Range should be 0-1
#'@param tol Numeric scalar. Tolerance level to add to biomass floor (Default = 1E-6)
#'@param plot A logical value specifying if the function should generate plots or
#'not. (Default = F).
#'
#'@importFrom magrittr %>%
#'
#'
#'@return Returns a data frame of species
#'
#'\item{species}{The common name of the species/functional group}
#'\item{Code}{Atlantis Code for species/functional group}
#'
#'@export
#'
#'
#'@examples
#'\dontrun{

#'}

diag_reasonability <- function(modelBiomass, initialYr=1964, speciesNames=NULL, realBiomass, nYrs = NULL,
                               bounds = c(1,1), tol = 1E-6, plot=F){

  # filter real biomass data
  realBiomass <- realBiomass %>%
    dplyr::select(YEAR,Species,value,Code,units) %>%
    dplyr::rename(year=YEAR,species=Species,code=Code)

  # filter model biomass and average over the year
  modelBiomass <- modelBiomass %>%
    dplyr::select(species, time, atoutput,code) %>%
    dplyr::rename(value=atoutput) %>%
    dplyr::mutate(yearStep= floor(time/365)) %>%
    dplyr::group_by(code,yearStep) %>%
    dplyr::summarise(aveBiomass = mean(value)) %>%
    dplyr::mutate(year = yearStep+initialYr,yearStep=NULL) %>%
    dplyr::ungroup() %>%
    tibble::as_tibble()

  # find final year of model run
  maxRuntime <- max(modelBiomass$year)

   # use default options
  if (is.null(speciesNames)) { # select all species
    speciesCodes <- unique(modelBiomass$code)
    speciesCodes <- speciesCodes[!is.na(speciesCodes)]
  } else { # find species code


  }

  if (is.null(nYrs)) { # use all time series
    filterTime <- initialYr
  } else { # last n years
    filterTime <- maxRuntime - nYrs + 1
  }

  # now find if model biomass falls withing real biomass * bounds
  # need to map model time to real time.
  # use %age of initial biomass if no realBiomass
  for (acode in speciesCodes[17]) {
    rb <- realBiomass %>%
      dplyr::filter(code == acode) %>%
      dplyr::filter(year > filterTime)

    if (nrow(rb)==0) { # no data

    } else {
      aSp <- rb %>% dplyr::distinct(species) %>% dplyr::pull()

      mb <- modelBiomass %>%
        dplyr::filter(code == acode) %>%
        dplyr::filter(year > filterTime)

    }
  }

   # mb %>% dplyr::left_join(.,rb)





}




