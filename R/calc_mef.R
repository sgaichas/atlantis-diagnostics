#' Calculate Modeling Efficiency (MEF)
#'
#' Internal function to calculate modeling efficiency
#'
#' @param observed observed values (data)
#' @param expected expected values (model)
#'
#'@return list of two elements
#'\item{mef}{model efficiency value: Range from -Inf to 1}
#'\item{n}{number of observations used to calculate mef}
#'
#'@noRd



calc_mef <- function(observed, expected) {

  #remove NAs. Exist only in observed data
  expected <- expected[!is.na(observed)]
  observed <- observed[!is.na(observed)]

  oo <- sum((observed-mean(observed))^2)
  po <- sum((expected - observed)^2)

  mef <- (oo-po)/oo

  # number of data point used in calculation
  n <- length(expected)

  return(list(mef=mef,n=n))
}
