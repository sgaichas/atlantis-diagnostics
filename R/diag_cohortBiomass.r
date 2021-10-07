#' Test for distribution of biomass over age class
#'
#' Determine which age classes contain the most biomass and compare to neighboring
#' age classes. A species fails the test is all biomass is contained in either
#' the smallest or largest age class
#'
#'@param mortFile A character string. Path to location of Mort.txt file.
#'@param agebiomind A character string. Path to location of AgeBiomIndx file.
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

diag_cohortBiomass <- function(mortFile,agebiomind,neusPriority) {

  cohortBiom <- read.csv(agebiomind,sep = " ", stringsAsFactors=FALSE, header=TRUE)
  mort <- read.csv(mortFile,sep = " ", stringsAsFactors=FALSE, header=TRUE)
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

  Code <- c("MAK","HER","WHK","BLF","WPF","SUF","WIF","WTF", "FOU", "HAL",	"PLA",	"FLA",	"BFT",	"TUN",	"BIL",	"MPF",	"BUT",	"BPF",	"ANC",	"GOO",	"MEN",	"FDE",	"COD",	"SHK",	"OHK",	"POL",	"RHK",	"BSB",	"SCU",	"TYL",	"RED",	"OPT",	"SAL",	"DRM",	"STB",	"TAU",	"WOL",	"SDF",	"FDF",	"HAD",	"YTF",	"DOG",	"SMO")
  numGroups <- length(Code)


  Max_Cohort <- c()
  Status <- c()
  Stability <- c()

  for (i in 1:numGroups) {
    groupName <- Code[i]
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


  diagnostics <- data.frame(Code,Status,Max_Cohort,Stability)
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
