#'Test for all functional group persistence
#'
#'@description \code{diag_persist} Averages over a year to see if mean B is 0,
#'assuming that initial (Time = 0) biomass does not count. Output is a list of
#'species that have crashed and an optional set of plots where 0 biomass is red.
#'
#'@details This should be run on an unfished, unperturbed run. Fishing or
#'perturbations may legitimately drive groups extinct. However, for our
#'historical period, we don’t expect anything to go extinct that is in the
#'system now, so this test includes historical fishing, physics, and
#'phytoplankton drivers.
#'
#'@param atBtxt A dataframe of total biomass of all groups over time, read in from
#'Atlantis ...BioInd.txt output using \code{atlantisom::load_bioind}.
#'@param fgs.names A vector of species names
#'@param plot A logical value specifying if the function should generate plots or
#'not. The default is \code{TRUE}.
#'
#'@importFrom magrittr %>%
#'@importFrom rlang .data
#'
#'@return Returns a dataframe of crashed species with columns "Failed" of species
#'names and "nYrs0B" with number of years having 0 biomass. Also plots all species
#'trajectories with red lines for crashed years unless plot=FALSE.
#'@export
#'
#'@author Sarah Gaichas
#'
#'@examples
#'\dontrun{
#'# Usage
#'diag_persist(testdiag$atBtxt, testdiag$fgs.names)
#'}
#'
diag_persist <- function(atBtxt, fgs.names, plot=TRUE){

  # need in annual units? Or fail when any output timestep below threshold?
  # make safe for migratory species, assume that over the course of the year mean B > 0.
  # assumes biomass never goes negative in atlantis

  crash <- atBtxt %>%
    dplyr::filter(.data$species %in% fgs.names) %>%
    dplyr::mutate(yr = ceiling(.data$time/365)) %>%
    dplyr::filter(.data$yr > 0) %>%
    dplyr::group_by(.data$species, .data$yr) %>%
    dplyr::summarise(meanB = mean(.data$atoutput)) %>%
    dplyr::filter(.data$meanB == 0)

  # flag any groups with any mean annual biomass below a threshold

  crashed <- crash %>%
    dplyr::group_by(.data$species) %>%
    dplyr::count()

  names(crashed) <- c("Failed", "nYrs0B")

  # visualize; hardcoded pages for ~89 group NEUS model

  if(plot){

    atBtxt$col <- cut(atBtxt$atoutput,
                      breaks = c(-Inf, 0, Inf),
                      labels = c("crashed", ">0 B"))

    plotB <-ggplot2::ggplot() +
      ggplot2::geom_line(data=atBtxt%>%dplyr::filter(.data$species %in% fgs.names),
                         ggplot2::aes(x=.data$time/365,y=.data$atoutput, color=.data$col),
                         alpha = 10/10) +
      ggthemes::theme_tufte() +
      ggplot2::theme(legend.position = "top") +
      ggplot2::labs(colour=g.name)

    print(plotB + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 1, scales="free"))
    print(plotB + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 2, scales="free"))
    print(plotB + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 3, scales="free"))
    print(plotB + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 4, scales="free"))
    print(plotB + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 5, scales="free"))
    print(plotB + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 6, scales="free"))
    print(plotB + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 7, scales="free"))
    print(plotB + ggforce::facet_wrap_paginate(~species, ncol=4, nrow = 3, page = 8, scales="free"))
  }

  return(as.data.frame(crashed))

}
