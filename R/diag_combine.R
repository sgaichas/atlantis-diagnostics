#' Combine all diagnostic output
#'
#' Using the output from all tests (persistence, stability, reasonability), they are combined
#' into a single data frame
#'
#' @param persistence Tibble. Output from \code{diag_persistence}
#' @param stability Tibble. Output from \code{diag_stability}
#' @param reasonability Tibble. Output from \code{diag_reasonability}
#'
#' @return A data frame
#' \item{species}{The common name of the species/functional group as described in Atlantis input file}
#' \item{code}{Atlantis Code for species/functional group}
#' \item{persistence}{Indicating pass or failure of persistence test}
#' \item{proportionInitBio}{Metric for persistence. The proportion of initial biomass at \code{tminimumBiomass}}
#' \item{stability}{Indicating pass or failure of stability test}
#' \item{relChange}{Metric for stability. Rate of increase relative to \code{t1Biomass}(\code{mtperyear}/\code{t1Biomass}). see \code{\url{diag_stability}} }
#' \item{reasonability}{Indicating pass or failure of reasonability test}
#' \item{}{}
#'
#'@export

diag_combine <- function(persistence, stability, reasonability){

  diagnostics <- persistence %>%
    dplyr::left_join(.,stability,by = c("code","species")) %>%
    dplyr::rename(persistence = pass.x, stability = pass.y) %>%
    dplyr::select(code, species, persistence, stability, proportionInitBio, relChange) %>%
    dplyr::left_join(.,reasonability,by=c("code","species")) %>%
    dplyr::rename(reasonability = pass) %>%
    dplyr::select(code, species, persistence,  proportionInitBio, stability, relChange,reasonability, maxExceedance) %>%
    dplyr::arrange(persistence, reasonability, stability)

  return(diagnostics)

}


