diag_persist <- function(atBtxt, fgs.names, plot=TRUE){
  
  # need in annual units? Or fail when any output timestep below threshold?
  # make safe for migratory species, assume that over the course of the year mean B > 0.
  # assumes biomass never goes negative in atlantis
  
  crash <- atBtxt %>%
    filter(species %in% fgs.names) %>%
    mutate(yr = ceiling(time/365)) %>%
    filter(yr > 0) %>%
    group_by(species, yr) %>%
    summarise(meanB = mean(atoutput)) %>%
    filter(meanB == 0)
  
  # flag any groups with any mean annual biomass below a threshold
  
  crashed <- crash %>%
    group_by(species) %>%
    count()
  
  names(crashed) <- c("Failed", "nYrs0B")
  
  # visualize; hardcoded pages for ~89 group NEUS model
  
  if(plot){
    
    atBtxt$col <- cut(atBtxt$atoutput,
                      breaks = c(-Inf, 0, Inf),
                      labels = c("crashed", ">0 B"))
    
    plotB <-ggplot() +
      geom_line(data=atBtxt%>%filter(species %in% fgs.names), 
                aes(x=time/365,y=atoutput, color=col),
                alpha = 10/10) +
      ggthemes::theme_tufte() +
      theme(legend.position = "top") +
      labs(colour=g.name)
    
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
