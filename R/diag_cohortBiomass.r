#' Test for distribution of biomass over age class
#'
#' Determine which age classes contain the most biomass and compare to neighboring
#' age classes. A species fails the test if all biomass is contained in either
#' the smallest or largest age class
#'
#'@param fgs A character string. Path to location of functional groups file.
#'@param mortality A character string. Path to location of Mort.txt file.
#'@param agebiomind A character string. Path to location of AgeBiomIndx file.
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for persistence.
#'(Default = NULL, uses all species with \code{IsTurnedOn=1} in \code{fgs} file)
#'@param neusPriority A character string. Path to location of Species Priorities file.
#'
#'@return
#'\item{code}{Atlantis species code}
#'\item{status}{Boolean indicating whether species passes test}
#'\item{maxCohort}{Age class with the highes biomass}
#'\item{stability}{}
#'\item{priority}{Scale defining priority of species in the model, High (H), Low(L) }
#'\item{fishing}{Scale defining level of fishing, 1 is highest}
#'
#' @export

diag_cohortBiomass <- function(fgs,
                               mortality,
                               agebiomind,
                               speciesCodes = NULL,
                               neusPriority) {


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


  cohortBiom <- utils::read.csv(agebiomind,sep = " ", stringsAsFactors=FALSE, header=TRUE)
  mort <- utils::read.csv(mortality,sep = " ", stringsAsFactors=FALSE, header=TRUE)
  mort <- dplyr::select(mort,dplyr::contains(".F"))
  mort_l20 <- dplyr::slice(mort,-c(1:34,55))
  meanMort <- dplyr::summarize_all(mort_l20,mean)


  neusPriority <- utils::read.csv(neusPriority,sep=",",stringsAsFactors=FALSE,header=TRUE) %>%
    dplyr::rename(priority = .data$priority.overall) %>%
    dplyr::select(.data$code,.data$priority)

  numRows <- nrow(cohortBiom)
  lastRow <- numRows - 1
  firstRow <- numRows - 100
  cohortBiom <- dplyr::slice(cohortBiom,firstRow:lastRow)

  numGroups <- length(speciesCodes)


  maxCohort <- c()
  pass <- c()
  stability <- c()

  for (i in 1:numGroups) {
    groupName <- speciesCodes[i]
    groupCohort <- dplyr::select(cohortBiom,contains(groupName))
    groupCohortMean <- dplyr::summarise_each(groupCohort, mean)
    maxCohortMean <- base::which.max(groupCohortMean)
    maxCohort <- c(maxCohort, maxCohortMean)
    if (maxCohortMean == 1 || maxCohortMean == 10) {
      pass <- c(pass,F)
    } else {
      pass <- c(pass, T)
    }
    maxMeanIndex <- base::which.max(groupCohortMean)
    maxMeanIndex <- maxMeanIndex[[1]]
    stabVal <- groupCohort[100,maxMeanIndex] / groupCohort[5,maxMeanIndex]

    if (stabVal < 0.75) {
      stability <- c(stability,"Declining")
    } else if (stabVal > 1.25) {
      stability <- c(stability, "Increasing")
    } else {
      stability <- c(stability, "Stable")
    }
  }


  diagnostics <- data.frame(code=speciesCodes,pass,maxCohort,stability)
  diagnostics <- dplyr::inner_join(diagnostics,neusPriority,by="code")
  diagnostics$fishing <- NULL

  numGroups <- nrow(diagnostics)
  for (i in 1:numGroups) {
    groupCol <- dplyr::select(meanMort,contains(diagnostics$code[i]))
    #print(groupCol[1,1])
    if (groupCol[1,1] >= 0.05) {
      diagnostics$fishing[i] <- 1
    } else if (groupCol[1,1] > 0.001) {
      diagnostics$fishing[i] <- 2
    } else if (groupCol[1,1] > 0.0001) {
      diagnostics$fishing[i] <- 3
    } else {
      diagnostics$fishing[i] <- 4
    }
  }

  diagnostics <- diagnostics %>%
    dplyr::arrange(.data$priority,.data$pass,.data$maxCohort,.data$stability,.data$fishing,.data$code) %>%
    tibble::as_tibble()

  return(diagnostics)
}
