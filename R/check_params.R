#' Check parameters exist
#'
#' Inspects biology parameter file to ensure all parameter values have been entered. Missing values will not
#' cause the model to crash but will produce unexpected results. Explain ...
#'
#' @param prmDir Character vector. Path to biology parameter file
#' @param fgs Data frame. Contents of the functional group file. \code{atlantisom::load_fgs}
#'
#'@importFrom magrittr %>%
#'
#'
#'@return Returns a data frame of species which do not meet defined persistence criteria.
#'
#'\item{a}{The common name of the species/functional group}
#'
#'@export
#'
#'
#'@examples
#'\dontrun{
#'}

check_params <- function(prmDir, fgs){

  a <- make.prm.attributes(prmDir, fgs) %>%
    tibble::aas_tibble()

  return(a)

}
