#'Test for functional group reasonability
#'
#' Inspects each time point. If at any point in time a species falls below or rises above predefined
#' bounds it is flagged. The term "reasonable" is based on survey data when available.
#'
#'@param fgs A character string. Path to location of functional groups file.
#'@param biomind A character string. Path to the BiomIndx.txt file.
#'@param initialYr Numeric Scalar. Year in which the model run was initiated. (Default = 1964)
#'@param startYr Numeric Scalar. Year in which the model finishes spin up period. (Default = 1998)
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for reasonableness
#'(Default = NULL, uses all species found in  \code{modelBiomass}. Species codes should be a subset of the Atlantis species codes
#'@param realBiomass A data frame. biomass time series (from assessments, stock SMART or otherwise) for species.
#' \code{realBiomass} should be in long format with column labels (YEAR,variable,value,Code,Species,isFishedSpecies)
#' variable should contain value = "biomass" (Biomass units should be in metric tonnes) and/or "var" if \code{useVariance = T}
#'@param useVariance Boolean. If to use variance estimates of biomass (included in \code{realBiomass}) (Default = F, reverts to using \code{surveyBounds} as upper and lover bounds)
#'@param nYrs Numeric scalar. Number of years from the end of the time series for which reasonableness is checked.
#' (Default = NULL, contemporary period from startYr is used)
#'@param surveyBounds Numeric vector. Size of 1x2 containing the values in which to multiple lower and upper bounds of observed data.
#'For example (Default = c(1,1)) indicating use of min and max of observed biomass
#'@param initBioBounds Numeric vector. Size of 1x2 containing lower and upper bound proportions used to scale initial biomass.
#'This is used for groups/species that dont have surveys
#'
#'
#'@return Returns a data frame of species
#'
#'\item{species}{The common name of the species/functional group as described in Atlantis input file}
#'\item{code}{Atlantis Code for species/functional group}
#'\item{initialBiomass}{Starting value of Biomass for species/functional group. From model output}
#'\item{minBiomass}{The smallest value of biomass observed in the run}
#'\item{maxBiomass}{The largest value of biomass observed in the run}
#'\item{propInitiBio}{The maxBiomass as a proportion of the inintial biomass}
#'\item{propBelowLower}{Magnitude of lower bound excedence as a proportion of lower threshold}
#'\item{propAboveUpper}{Magnitude of upper bound excedence as a proportion of upper threshold}
#'\item{maxExceedance}{max of \code{propBelowLower} and \code{propAboveUpper}}
#'\item{t1}{The first year reasonableness was not met}
#'\item{tn}{The last year reasonableness was not met}
#'\item{nts}{The total number of years that reasonableness was not met}
#'\item{modelSkill}{Measure of the goodness of fit}
#'\item{nObs}{number of observations used in the goodness of fit}
#'\item{pass}{Boolean. Indicates if the Atlantis group passed reasonability}
#'\item{test}{Character. Indicate how \code{pass} was determined. If the Atlantis species is fished and the species is a stand alone species in the model
#'then survey data is used to assess reasonability. Otherwise deviations from initial biomass is used}
#'
#'
#'@family diags
#'
#'@export
#'
#'@importFrom magrittr %>%
#'@importFrom rlang .data
#'
#'@examples
#'\dontrun{
#'# Declare paths to files required
#'
#' biomind <- paste("Full path to file","xxxBiomIndx.txt")
#' fgs <- paste("Full path to file","functionalGroups.csv")
#'
#' # read in survey biomass and convert to metric tons
#' realBiomass <- readRDS(paste0(dataDir,"sweptAreaBiomassNEUS.rds")) %>%
#'       dplyr::filter(variable %in% c("tot.biomass")) %>%
#'       dplyr::mutate(value=ifelse(grepl("kg$",units),value/1000,value)) %>%
#'       dplyr::select(-units)
#'
#' # Perform reasonability test on all species/groups using the last 20 years of the run.
#' # Allow species with data to be bounded by 100 x max observed biomass and 1 x min observed biomass.
#' # For species without data allow model biomass to lie between 0.5 and 2 times initial biomass
#'
#' diag_reasonability(fgs, biomind, initialYr = 1964, realBiomass=realBiomass,
#' surveyBounds = c(1,100), initBioBounds = c(0.5,2))
#'
#' # Only perform test on herring and white hake.
#' diag_reasonability(fgs, biomind, initialYr = 1964, speciesCodes =c("MAK","WHK"),
#'  realBiomass=realBiomass, surveyBounds = c(1,100), initBioBounds = c(0.5,2))
#'}

diag_reasonability <- function(fgs,
                               biomind,
                               initialYr=1964,
                               startYr = 1998,
                               speciesCodes=NULL,
                               realBiomass,
                               useVariance=F,
                               nYrs = NULL,
                               surveyBounds = c(1,1),
                               initBioBounds = c(0.5,10)){

  # Check realBiomass for variance values
  if (useVariance == T) {
    if (!(any(realBiomass %>% dplyr::distinct(variable) == "var"))) {
      stop(" To use variance estimates in reasonability bounds the values of
           variance need to be included in the realBiomass argument")
    }

  }

  ################################################
  ########### model output #######################
  ################################################
  # read in biomass data and qualify species Codes
  biom <- get_model_biomass(fgs,biomind,speciesCodes)
  modelBiomass <- biom$modelBiomass
  speciesCodes <- biom$speciesCodes

  # list of Atlantis species and Atlantis codes. Pulled from model output
  species <- modelBiomass %>%
    dplyr::select(.data$species,.data$code) %>%
    dplyr::filter(!is.na(.data$code)) %>%
    dplyr::distinct(.data$species,.data$code)

  # pull initial biomass for each species/group
  initialBiomass <- modelBiomass %>%
    dplyr::group_by(.data$species) %>%
    dplyr::mutate(initialBiomass = dplyr::first(.data$atoutput)) %>%
    dplyr::select(.data$species,.data$code,.data$initialBiomass) %>%
    dplyr::distinct(.data$species,.data$code,.data$initialBiomass)

  # filter model biomass and average over the year
  # join with initial biomass
  modelBiomass <- modelBiomass %>%
    dplyr::select(.data$time, .data$atoutput, .data$code, .data$species) %>%
    dplyr::rename(value=.data$atoutput) %>%
    dplyr::mutate(yearStep= floor(.data$time/365)) %>%
    dplyr::group_by(.data$code,.data$yearStep,.data$species) %>%
    dplyr::summarise(modelBiomass = mean(.data$value),.groups="drop") %>%
    dplyr::mutate(year = .data$yearStep+initialYr,yearStep=NULL) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(.,initialBiomass,by=c("code","species"))

  ################################################
  ############## observed data ###################
  ################################################

  # filter real biomass data to select only species that are fished in the atlantis model
  # then remove sharks, whales, birds etc if present.

  speciesBiomass <- realBiomass %>%
    dplyr::filter(.data$isFishedSpecies == T) %>%
    dplyr::filter(!grepl("shark",.data$Species)) %>%
    dplyr::filter(!grepl("whale",.data$Species)) %>%
    dplyr::filter(!grepl("bird",.data$Species)) %>%
    dplyr::filter(!grepl("salmon",.data$Species)) %>%
    dplyr::select(.data$YEAR,.data$variable,.data$value,.data$Code,.data$Species) %>%
    dplyr::rename(year=.data$YEAR,code=.data$Code)

  #reorder observed data and make wide format
  observedData <- speciesBiomass %>%
    dplyr::group_by(.data$year,.data$variable,.data$code) %>%
    dplyr::summarise(value = sum(value),.groups = "drop") %>% # aggregate over species with same functional group
    tidyr::pivot_wider(.,names_from = .data$variable,values_from = .data$value)


  # find final year of model run
  maxRuntime <- max(modelBiomass$year)

  # determine time frame in which to perform "test"
  if (is.null(nYrs)) { # use all time series
    filterTime <- startYr #initialYr
  } else { # last n years
    filterTime <- max(startYr,maxRuntime - nYrs + 1)
  }

  # now find if model biomass falls within real biomass * bounds
  # need to map model time to real time.
  # use %age of initial biomass if no real Biomass
  reasonable <- NULL
  # loop over species
  for (acode in speciesCodes) {

    # filter real/observed data by species and the time frame
    rb <- observedData %>%
      dplyr::filter(.data$code == acode) %>%
      dplyr::filter(.data$year >= filterTime)


    # if data available for species code then get model data for time range
    if (nrow(rb)==0) { # no data
      # take %age of initialbiomass

      # compare model biomass output against scaled initial biomass
      # calculate goodness of fit (mef)
      mb <- modelBiomass %>%
        dplyr::filter(.data$code == acode) %>%
        dplyr::filter(.data$year >= filterTime) %>%
        dplyr::mutate(reasonable = (.data$modelBiomass < .data$initialBiomass*initBioBounds[2]) &
                        (.data$modelBiomass >= .data$initialBiomass*initBioBounds[1])) %>%
        dplyr::mutate(modelSkill = calc_mef(.data$modelBiomass,.data$initialBiomass)$mef) %>%
        dplyr::mutate(minBiomass = min(.data$modelBiomass,na.rm=T)) %>%
        dplyr::mutate(maxBiomass = max(.data$modelBiomass,na.rm=T)) %>%
        dplyr::mutate(propInitBio = .data$maxBiomass/.data$initialBiomass) %>%
        dplyr::mutate(propAboveUpper = .data$maxBiomass/(.data$initialBiomass*initBioBounds[2]) - 1) %>%
        dplyr::mutate(propBelowLower = -.data$minBiomass/(.data$initialBiomass*initBioBounds[1]) + 1) %>%
        dplyr::group_by(.data$code) %>%
        dplyr::mutate(maxExceedance = max(.data$propBelowLower,.data$propAboveUpper)) %>%
        dplyr::ungroup()

      test <- "initialBio"

    } else { # if complete data present

      if (!useVariance) { # use surveybounds argument to determine reasonability
        # Upper bound is maximum observed biomass (over time series) * surveyBounds[2]
        # compare model biomass against scaled survey estimates
        # calculate goodness of fit (mef)
        mb <- modelBiomass %>%
          dplyr::filter(.data$code == acode) %>%
          dplyr::filter(.data$year >= filterTime) %>%
          dplyr::left_join(.,rb,by=c("code"="code","year"="year")) %>%
          dplyr::mutate(reasonable = (.data$modelBiomass < max(.data$biomass,na.rm=T)*surveyBounds[2]) &
                          (.data$modelBiomass > min(.data$biomass,na.rm=T)*surveyBounds[1])) %>%
          dplyr::mutate(modelSkill = calc_mef(.data$biomass,.data$modelBiomass)$mef) %>%
          dplyr::mutate(minBiomass = min(.data$modelBiomass)) %>%
          dplyr::mutate(maxBiomass = max(.data$modelBiomass)) %>%
          dplyr::mutate(propInitBio = .data$maxBiomass/.data$initialBiomass) %>%
          dplyr::mutate(propAboveUpper = .data$maxBiomass/(max(.data$biomass,na.rm=T)*surveyBounds[2]) - 1) %>%
          dplyr::mutate(propBelowLower = -.data$minBiomass/(min(.data$biomass,na.rm=T)*surveyBounds[1]) + 1) %>%
          dplyr::group_by(.data$code) %>%
          dplyr::mutate(maxExceedance = max(.data$propBelowLower,.data$propAboveUpper)) %>%
          dplyr::ungroup()
      } else { # bounds based on variance of estimate of biomass
        # Upper bound is maximum (observed biomass * 3*max sd observed) over time series
        mb <- modelBiomass %>%
          dplyr::filter(.data$code == acode) %>%
          dplyr::filter(.data$year >= filterTime) %>%
          dplyr::left_join(.,rb,by=c("code"="code","year"="year")) %>%
          dplyr::mutate(upperBound = .data$biomass + 3* sqrt(.data$var)) %>%
          dplyr::mutate(lowerBound = .data$biomass - 3* sqrt(.data$var)) %>%
          dplyr::mutate(reasonable = (.data$modelBiomass < max(.data$upperBound,na.rm=T)) &
                          (.data$modelBiomass > max(0,.data$lowerBound,na.rm=T))) %>%
          dplyr::mutate(modelSkill = calc_mef(.data$biomass,.data$modelBiomass)$mef) %>%
          dplyr::mutate(minBiomass = min(.data$modelBiomass)) %>%
          dplyr::mutate(maxBiomass = max(.data$modelBiomass)) %>%
          dplyr::mutate(propInitBio = .data$maxBiomass/.data$initialBiomass) %>%
          dplyr::mutate(propAboveUpper = .data$maxBiomass/(max(.data$upperBound,na.rm=T)) - 1) %>%
          dplyr::mutate(propBelowLower = -.data$minBiomass/(max(0,.data$lowerBound,na.rm=T)) + 1) %>%
          dplyr::group_by(.data$code) %>%
          dplyr::mutate(maxExceedance = max(.data$propBelowLower,.data$propAboveUpper)) %>%
          dplyr::ungroup()

      }

      test <- "data"
    }

    # filter all years failing test
    out <- mb %>%
      dplyr::filter(.data$reasonable == F)


    if(nrow(out) == 0) { # all pass
      # create output
      out <- mb %>%
        dplyr::select(.data$code,.data$species,.data$initialBiomass,.data$modelSkill,.data$minBiomass,.data$maxBiomass,
                      .data$propInitBio,.data$propBelowLower,.data$propAboveUpper,.data$maxExceedance) %>%
        dplyr::distinct() %>%
        dplyr::mutate(t1 = NA) %>%
        dplyr::mutate(tn = NA) %>%
        dplyr::mutate(nts = 0) %>%
        dplyr::mutate(pass = T) %>%
        dplyr::mutate(test = test)
    } else { # some failure
      out <- mb %>%
        dplyr::filter(.data$reasonable == F) %>%
        dplyr::mutate(t1 = min(.data$year)) %>%
        dplyr::mutate(tn = max(.data$year)) %>%
        dplyr::group_by(.data$code,.data$species,.data$initialBiomass,.data$t1,.data$tn,.data$modelSkill,.data$minBiomass,.data$maxBiomass,
                        .data$propInitBio,.data$propBelowLower,.data$propAboveUpper,.data$maxExceedance) %>%
        dplyr::summarize(nts = sum(!.data$reasonable),.groups="drop") %>%
        dplyr::mutate(pass = dplyr::if_else(.data$nts>0,F,T)) %>%
        dplyr::mutate(test = test)

    }

    # bind to main output data frame
    reasonable <- rbind(reasonable,out)

  }
  reasonable <- reasonable %>%
    dplyr::arrange(.data$pass,dplyr::desc(.data$maxExceedance))

  return(reasonable)


}


