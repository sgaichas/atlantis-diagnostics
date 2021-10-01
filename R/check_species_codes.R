#' Check for valid species codes
#'
#' Uses complete list of atlantis codes to ensure codes passed to test are part of the list
#'
#' @param modelBiomass Data frame. Output from \code{atlantisom::load_bioind}
#' @param speciesCodes Character vector. Vector of species Codes to verify against complete list
#'
#' @return Vector of valid Codes
#'
#' @importFrom rlang .data
#'
#' @noRd

check_species_codes <- function(modelBiomass,speciesCodes){

  # use default options
  atlantisCodes <- modelBiomass %>%
    dplyr::filter(!is.na(.data$code)) %>%
    dplyr::distinct(.data$code) %>%
    dplyr::pull()
  # use default options
  if (is.null(speciesCodes)) { # select all species
    speciesCodes <- atlantisCodes
  } else {
    # remove NA's check for codes not in atlantis
    speciesCodes <- speciesCodes[!is.na(speciesCodes)]
    # check to make sure no non atlantis codes
    invalidCodes <- base::setdiff(speciesCodes,atlantisCodes)

    if (!(length(invalidCodes)==0)){
      stop("Invalid Atlantis group codes: ",paste0(invalidCodes,collapse=", "))
    }
  }

  return(speciesCodes)

}
