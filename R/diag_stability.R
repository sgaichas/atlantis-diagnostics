#'Test for functional group stability
#'
#'\code{diag_stability} determines whether the model reaches a steady state
#'over the last n years of a run. Stability is loosely defined as a species/groups biomass reaching a
#'stable level (measured as having no trend in a defined time range).
#'
#'
#'@param modelBiomass A data frame. Total biomass of all groups over time, read in from
#'Atlantis ...BioInd.txt output using \code{atlantisom::load_bioind}.
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for stability.
#'(Default = NULL, uses all species found in  \code{modelBiomass})
#'@param nYrs Numeric scalar. Number of years from the end of the time series that stability must occur.
#' (Default = 20 years)
#'@param sigTest Numeric Scalar. alpha level to use to test for slope significance (Default = 0.05)
#'@param plot Logical. Specifying whether the function should generate plots or not. (Default is F).
#'
#'@importFrom magrittr %>%
#'
#'@return Returns a data frame of all species and how they measure up against the stability criterion
#'\item{code}{Atlantis Code for species/functional group}
#'\item{species}{The common name of the species/functional group}
#'\item{slope}{Double. The value of the slope parameter (time) }
#'\item{pValue}{Double. The p-value associated with the slope parameter (time)}
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
#' # Perform stability test on all species/groups using the last 20 years of the run
#' diag_stability(modelBiomass, nYrs = 20)
#'
#' # Only perform test on herring and white hake.
#' # Require stability over the last 10 years of the run and and use alpha = 0.1 as tests significance level
#' diag_stability(modelBiomass, speciesCodes=c("HER","WHK"), nYrs = 10, sigTest = 0.1)
#'}

diag_stability <- function(modelBiomass, speciesCodes, nYrs = 20, sigTest = 0.05, plot=F){

  # find final time step value
  maxRuntime <- max(modelBiomass$time)

  # check for valid species Codes & clean
  speciesCodes <- check_species_codes(modelBiomass,speciesCodes)


  if (is.null(nYrs)) { # use all time series
    filterTime <- 0
  } else { # last n years
    filterTime <- maxRuntime - (365*nYrs)
  }

  # select species, timframe, group, nest
  stable <- modelBiomass %>%
    dplyr::filter(time >= filterTime) %>%
    dplyr::filter(code %in% speciesCodes) %>%
    dplyr::select(code, species,time,atoutput) %>%
    dplyr::group_by(code,species) %>%
    tidyr::nest()

  # fit linear model (biomass = alpha + time.beta) to each species over the time frame specified
  # then extract slope and significance from model fit
  stability <- stable %>%
    dplyr::mutate(model = purrr::map(data,fitlm)) %>%
    dplyr::transmute(species,slope = purrr::map_dbl(model,coefs),pValue=purrr::map_dbl(model,pVals)) %>%
    dplyr::mutate(pass = pValue > sigTest) %>%
    dplyr::ungroup()


  stability <- stability %>%
    dplyr::arrange(pass,pValue)

  return(stability)



}


## Helper functions

coefs <- function(model){
  coefficients(model)[["time"]]
}

pVals <- function(model){
  summary(model)$coefficients[,4][["time"]]
}

fitlm <- function(df){
  lm(atoutput ~ time, data=df)
}







