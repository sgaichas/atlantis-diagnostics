#' Helper functions
#'
#'
#'@noRd

coefs <- function(model){
  stats::coefficients(model)[["year"]]
}

#'@noRd
pVals <- function(model){
  base::summary(model)$coefficients[,4][["year"]]
}

#'@noRd
fitlm <- function(df){
  stats::lm(meanBio ~ year, data=df)
}

#'@noRd
fittedVal <- function(model) {
  # first fitted value
  stats::fitted.values(model)[[1]]
}

#'@noRd
meanData <- function(df){
  base::mean(df$meanBio)
}
