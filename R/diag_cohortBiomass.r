#' Test for distribution of biomass over age class
#'
#' Determine which age classes contain the most biomass and compare to neighboring
#' age classes. A species fails the test is all biomass is contained in either
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
#'\item{Code}{Atlantis species code}
#'\item{Status}{Boolean indicating whether species passes test}
#'\item{Max_Cohort}{Age class with the highes biomass}
#'\item{Stability}{}
#'\item{Priority}{Scale defining priority of species in the model}
#'\item{Fishing}{Scale defining level of fishing, 1 is highest}
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


  cohortBiom <- read.csv(agebiomind,sep = " ", stringsAsFactors=FALSE, header=TRUE)
  mort <- read.csv(mortality,sep = " ", stringsAsFactors=FALSE, header=TRUE)
  mort <- dplyr::select(mort,dplyr::contains(".F"))
  mort_l20 <- dplyr::slice(mort,-c(1:34,55))
  meanMort <- dplyr::summarize_all(mort_l20,mean)

  #setwd("C:/Users/robert.gamble/Desktop/Atlantis_1_5/neus-atlantis/diagnostics")

  neusPriority <- read.csv(neusPriority,sep=",",stringsAsFactors=FALSE,header=TRUE) %>%
    dplyr::rename(Code = .data$code, Priority = .data$priority.overall)
  neusPriority <- dplyr::select(neusPriority, c(Code,Priority))
  numRows <- nrow(cohortBiom)
  lastRow <- numRows - 1
  firstRow <- numRows - 100
  cohortBiom <- dplyr::slice(cohortBiom,firstRow:lastRow)

  #Code <- c("MAK","HER","WHK","BLF","WPF","SUF","WIF","WTF", "FOU", "HAL",	"PLA",	"FLA",	"BFT",	"TUN",	"BIL",	"MPF",	"BUT",	"BPF",	"ANC",	"GOO",	"MEN",	"FDE",	"COD",	"SHK",	"OHK",	"POL",	"RHK",	"BSB",	"SCU",	"TYL",	"RED",	"OPT",	"SAL",	"DRM",	"STB",	"TAU",	"WOL",	"SDF",	"FDF",	"HAD",	"YTF",	"DOG",	"SMO")
  numGroups <- length(speciesCodes)


  Max_Cohort <- c()
  Status <- c()
  Stability <- c()

  for (i in 1:numGroups) {
    groupName <- speciesCodes[i]
    groupCohort <- dplyr::select(cohortBiom,contains(groupName))
    groupCohortMean <- dplyr::summarise_each(groupCohort, mean)
    maxCohortMean <- base::which.max(groupCohortMean)
    Max_Cohort <- c(Max_Cohort, maxCohortMean)
    if (maxCohortMean == 1 || maxCohortMean == 10) {
      Status <- c(Status,"FAIL")
    } else {
      Status <- c(Status, "PASS")
    }
    maxMeanIndex <- base::which.max(groupCohortMean)
    maxMeanIndex <- maxMeanIndex[[1]]
    stabVal <- groupCohort[100,maxMeanIndex] / groupCohort[5,maxMeanIndex]

    if (stabVal < 0.75) {
      Stability <- c(Stability,"  Declining")
    } else if (stabVal > 1.25) {
      Stability <- c(Stability, "  Increasing")
    } else {
      Stability <- c(Stability, "  Stable")
    }
  }


  diagnostics <- data.frame(Code=speciesCodes,Status,Max_Cohort,Stability)
  diagnostics <- dplyr::inner_join(diagnostics,neusPriority,by="Code")
  diagnostics$Fishing <- "4 - None"

  numGroups <- nrow(diagnostics)
  for (i in 1:numGroups) {
    groupCol <- dplyr::select(meanMort,contains(diagnostics$Code[i]))
    #print(groupCol[1,1])
    if (groupCol[1,1] >= 0.05) {
      diagnostics$Fishing[i] <- "1"
    } else if (groupCol[1,1] > 0.001) {
      diagnostics$Fishing[i] <- "2"
    } else if (groupCol[1,1] > 0.0001) {
      diagnostics$Fishing[i] <- "3"
    } else {
      diagnostics$Fishing[i] <- "4"
    }
  }

  diagnostics <- dplyr::arrange(diagnostics,Priority,Status,Max_Cohort,Stability,Fishing,Code)
  rownames(diagnostics) <-c()
  return(diagnostics)
}
