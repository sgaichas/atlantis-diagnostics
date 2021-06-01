#'Test for functional group stability
#'
#'\code{diag_stability} determines whether the model reaches a steady state
#'over the last n years of a run. Stability is loosely defined as a species/groups biomass reaching a
#'stable level (measured as having tolerable trend in a defined time range).
#'
#'
#'@param modelBiomass A data frame. Total biomass of all groups over time, read in from
#'Atlantis ...BioInd.txt output using \code{atlantisom::load_bioind}.
#'@param initialYr Numeric Scalar. Year in which the model run was initiated. (Default = 1964)
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for stability.
#'(Default = NULL, uses all species found in  \code{modelBiomass})
#'@param nYrs Numeric scalar. Number of years from the end of the time series that stability must occur.
#' (Default = 20 years)
#'@param relChangeThreshold Numeric Scalar. Maximum magnitude of relative change of slope (Default = 0.05)
#'@param plot Logical. Specifying whether the function should generate plots or not. (Default is F).
#'
#'@importFrom magrittr %>%
#'
#'@return Returns a data frame of all species and how they measure up against the stability criterion
#'\item{code}{Atlantis Code for species/functional group}
#'\item{species}{The common name of the species/functional group}
#'\item{t1Fit}{Value of fitted biomass for the first of year data used in the fit}
#'\item{mtPerYear}{Double. The value of the slope parameter (year) }
#'\item{relChange}{Rate of increase relative to \code{t1Fit}(\code{mtPerYear}/\code{t1Fit})}
#'\item{aveBio}{mean biomass for the last \code{nYrs} years }
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
#' Note: annual biomass is used in fitting. Species with mean annual biomass < 1 metric ton over the last
#' n years of the run are not considered stable. They are reported to Fail the test and NaNs returned
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
#' diag_stability(modelBiomass, speciesCodes=c("HER","WHK"), nYrs = 10, relChangeThreshold = 0.01)
#'}

diag_stability <- function(modelBiomass, initialYr = 1964, speciesCodes, nYrs = 20, relChangeThreshold = 0.01, plot=F){

  # check for valid species Codes & clean
  speciesCodes <- check_species_codes(modelBiomass,speciesCodes)

  # determine time frame in which to perform "test"
  # filter model biomass and average over the year
  # join with initial biomass
  modelBiomass <- modelBiomass %>%
    dplyr::select(time, atoutput, code, species) %>%
    dplyr::rename(value=atoutput) %>%
    dplyr::mutate(yearStep= floor(time/365)) %>%
    dplyr::group_by(code,yearStep,species) %>%
    dplyr::summarise(meanBio = mean(value),.groups="drop") %>% # average over year, removes seasonal cycle
    dplyr::mutate(year = yearStep+initialYr,yearStep=NULL) %>%
    dplyr::ungroup()

  if (is.null(nYrs)) { # use all time series
    filterTime <- initialYr
  } else { # last n years
    filterTime <- max(modelBiomass$year) - nYrs
  }

  # select species, time frame, group, nest
  stable <- modelBiomass %>%
    dplyr::filter(year > filterTime) %>%
    dplyr::filter(code %in% speciesCodes) %>%
    dplyr::group_by(code,species) %>%
    #dplyr::mutate(t1Biomass = dplyr::first(meanBio)) %>%
    dplyr::group_by(code,species) %>%
    tidyr::nest() # produces a column called data

  # fit linear model (biomass = alpha + year.beta) to each species over the time frame specified
  # then extract slope and significance from model fit

  stability <- stable %>%
    dplyr::mutate(model = purrr::map(data,fitlm)) %>%
    #dplyr::transmute(species,mtPerYear = purrr::map_dbl(model,coefs),pValue=purrr::map_dbl(model,pVals)) %>%
    dplyr::transmute(species,t1Fit = purrr::map_dbl(model,fittedVal),
                     mtPerYear = purrr::map_dbl(model,coefs),
                     aveBio = purrr::map_dbl(data,meanData)) %>%
    #dplyr::mutate(pass = pValue > sigTest) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(relChange = mtPerYear/t1Fit) %>%
    dplyr::mutate(sign = base::sign(mtPerYear)) %>%
    dplyr::mutate(pass = abs(relChange) < relChangeThreshold) %>%
    # average biomass < 1 metric ton
    dplyr::mutate(pass = dplyr::if_else(aveBio < 1,F,pass)) %>%
    dplyr::mutate(sign = dplyr::if_else(aveBio < 1,NaN,sign)) %>%
    dplyr::mutate(relChange = dplyr::if_else(aveBio < 1,NaN,relChange)) %>%
    dplyr::mutate(t1Fit = dplyr::if_else(aveBio < 1,NaN,t1Fit)) %>%
    dplyr::mutate(mtPerYear = dplyr::if_else(aveBio < 1,NaN,mtPerYear))


  stability <- stability %>%
    dplyr::arrange(pass,desc(abs(relChange)))

  return(stability)



}


## Helper functions

coefs <- function(model){
  coefficients(model)[["year"]]
}

pVals <- function(model){
  summary(model)$coefficients[,4][["year"]]
}

fitlm <- function(df){
  lm(meanBio ~ year, data=df)
}

fittedVal <- function(model) {
  # first fitted value
  fitted.values(model)[[1]]
}

meanData <- function(df){
  mean(df$meanBio)
}






