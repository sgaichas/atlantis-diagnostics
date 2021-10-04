#' Chunk several diag_functions share
#'
#' Reads in biomind data and checks that user species species codes are species
#' codes used in the model
#'
#'@param fgs A character string. Path to location of functional groups file.
#'@param biomind A character string. Path to the BiomIndx.txt file.
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for persistence.
#'
#' @noRd
#'

get_model_biomass <- function(fgs,biomind,speciesCodes){

  # get species codes
  allCodes <- atlantistools::get_turnedon_acronyms(fgs)
  if(!is.null(speciesCodes)) { # user supplied codes
    # check to see if codes are valid model codes
    invalidCodes <- base::setdiff(speciesCodes,allCodes)
    if (!(length(invalidCodes)==0)){
      stop("Invalid Atlantis group codes: ",paste0(invalidCodes,collapse=", "))
    }
  } else { # use all codes
    speciesCodes <- allCodes
  }

  # Need common name of species for output. read in fgs file and select common name and code
  functionalGps <- atlantistools::load_fgs(fgs)
  speciesNames <- functionalGps %>%
    dplyr::select(.data$Code,.data$Name)

  # get biomass from biomIndx file
  modelBiomass <- atlantistools::load_txt(biomind)  %>%
    dplyr::filter(.data$code %in% speciesCodes) %>%
    dplyr::left_join(., speciesNames, by=c("code"="Code")) %>%
    dplyr::rename(species = .data$Name) %>%
    dplyr::relocate(.data$species, .before = .data$time) %>%
    dplyr::arrange(.data$species,.data$time)

  return(list(modelBiomass = modelBiomass,speciesCodes = speciesCodes))

}
