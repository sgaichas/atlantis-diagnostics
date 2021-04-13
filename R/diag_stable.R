#'Test for model stability
#'
#'@description \code{diag_stable} determines whether the model reaches steady state
#'for the last ~20 years of an unperturbed, unfished ~100 year run.
#'
#'@details This test evaluates the last 20 years of the input run for a significant
#'slope using the \code{geom_gls} function from the package \code{ecodata} at
#'\link{https://github.com/NOAA-EDAB/ecodata}. The function `diag_stable` also
#'checks for persistence and flags groups that have gone to 0 biomass--they don't
#'count as stable. It takes the same inputs as `diag_persist` with the addition of
#'the \code{runpar} object, and defaults to plotting, which can be turned off with
#'`plot=FALSE` in the function call.
#'
#'@param atBtxt A dataframe of total biomass of all groups over time, read in from
#'Atlantis ...BioInd.txt output using \code{atlantisom::load_bioind}.
#'@param fgs.names A vector of species names
#'@param runpar A list of Atlantis run parameters read in from the .xml file using
#'\code{atlantisom::load_runprm}
#'@param plot A logical value specifying if the function should generate plots or
#'not. The default is \code{TRUE}.
#'
#'@importFrom magrittr %>%
#'
#'@return Returns a dataframe of crashed species with columns "Failed" of species
#'names and "nYrs0B" with number of years having 0 biomass. Also plots all species
#'trajectories with orange (positive) or purple (negative) trend lines for unstable
#'groups over the last 20 years unless plot=FALSE.
#'@noRd
#'
#'@author Sarah Gaichas
#'
#'@examples
#'\dontrun{
#'# Usage
#'diag_stable(testdiag$atBtxt, testdiag$fgs.names, testdiag$runpar)
#'}

#'
diag_stable <- function(atBtxt, fgs.names, runpar, plot=TRUE){

  # look for non-significant slope over last 20 years? 30 years would be better
  nlast <- 20

  startlast <- floor(runpar$nyears)-nlast

  #warning, was getting different behavior with this code on my linux machine

  stable <- atBtxt%>%filter(species %in% fgs.names) %>%
    #filter(time %in% seq(startlast*365, floor(runpar$nyears)*365, by=365)) %>%
    filter(time %in% seq(startlast*365, floor(runpar$nyears)*365, by=runpar$outputstep)) %>%
    group_by(species) %>%
    ggplot(aes(x=time/365, y=atoutput)) +
    ggthemes::theme_tufte() +
    geom_line(data=atBtxt%>%filter(species %in% fgs.names),
              aes(x=time/365, y=atoutput),
              alpha = 5/10) +
    geom_gls(warn = FALSE)


  print(stable + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 1, scales="free"))
  print(stable + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 2, scales="free"))
  print(stable + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 3, scales="free"))
  print(stable + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 4, scales="free"))
  print(stable + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 5, scales="free"))
  print(stable + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 6, scales="free"))
  print(stable + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 7, scales="free"))
  print(stable + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 8, scales="free"))

  # what I need is to extract the output of geom_gls that is NULL (no significant trend)

  # this may work eventually but blows up R right now, fix later
  # # get results from ecotrend package
  # stabletest <- atBtxt %>%
  #   filter(species %in% fgs.names) %>%
  #   #mutate(yr = ceiling(time/365)) %>%
  #   filter(time %in% seq(startlast*365, floor(runpar$nyears)*365, by=365)) %>%
  #   group_by(species) %>%
  #   dplyr::do(ecotrend::glsMs(data = ., formula = time ~ atoutput))
  #
  #ecotrend::glsMs(data = ., formula = yr ~ biomass)

  #ecotrend::glsMs(yr ~ biomass, stabletest)

  #test this run for persistence too

  # flag any groups with any mean annual biomass below a threshold


}
