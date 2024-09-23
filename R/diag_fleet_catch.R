#'Test fleet catch
#'
#'\code{diag_fleet_catch} determines whether the simulated fleets catcht the right magnitude and spatial distribution of catch
#'over the last n years of a run.
#'
#'
#'@param fgs A character string. Path to location of functional groups file.
#'@param fishery.prm A character string. Path to the fishery definitions file.
#'@param catch.file A character string. Path to catch nc file
#'@param bgm A character string. Path to the bgm file.
#'@param catch.ref A data.frame containing reference catch by fleet data (species|fleet|polygon|ref.value)
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for stability.
#'@param nYrs Numeric scalar. Number of years from the end of the time series that stability must occur.
#'@param relChangeThreshold Numeric Scalar. Maximum magnitude of relative change of slope (Default = 0.01)
#'@param min.dist Numeric Scalar. Maximum distance between model and reference center of gravity
#'
#'@return Returns a data frame of all species and how they measure up against the catch by fleet criteria
#'\item{species}{The common name of the species/functional group}
#'\item{fleet}{The fleet name from the fisheries file}
#'\item{catch.model}{Double. The mean catch in mT over the last nYrs of the model }
#'\item{catch.ref}{Double. The reference catch (mT)}
#'\item{cog.x.model}{Double. The mean x-coordinate center of gravity from the model}
#'\item{cog.x.ref}{Double. The mean x-coordinate center of gravity from the reference}
#'\item{cog.y.model}{Double. The mean y-coordinate center of gravity from the model}
#'\item{cog.y.ref}{Double. The mean y-coordinate center of gravity from the reference}
#'\item{dist}{Double. The mean distance between the center of gravity between the model and reference}
#'\item{catch.magnitude}{Logical. Is the relative difference between model and refernece catch within +/- relChangeThreshold?}
#'\item{catch.dist}{Logical. Is the distance between center of gravity betwene model and reference less than min.dist?}
#'
#'@family diags
#'
#'@export
#'
#'@importFrom magrittr %>%

diag_fleet_catch <- function(fgs,
                             fishery.prm,
                             catch.file,
                             bgm,
                             catch.ref,
                             speciesCodes = NULL,
                             nYrs = 20,
                             min.dist = 100,
                             relChangeThreshold = 0.01){

  boxes = atlantistools::convert_bgm(bgm)%>%
    dplyr::distinct(polygon,inside_lat,inside_long)

  catch.fleet =atlantisprocessing::process_catch_fleet(fishery.prm = fishery.prm,
                                          catch = catch.file,
                                          groups.file = fgs)

  if(!is.null(speciesCodes)){
    fgs.df = read.csv(fgs,as.is =T)

    spp.match  = fgs.df$LongName[which(fgs.df$Code %in% speciesCodes)]

    catch.fleet = catch.fleet %>%
      dplyr::filter(species %in% spp.match)

    catch.ref = catch.ref %>%
      dplyr::filter(species %in% spp.match)
  }
  max.yr = max(catch.fleet$time)

  #mnagnitude
  catch.mag.model =catch.fleet %>%
    dplyr::filter(time >= (max.yr - nYrs)) %>%
    dplyr::group_by(species,fleet,time)%>%
    dplyr::summarise(catch.model = sum(atoutput,na.rm=T))%>%
    dplyr::group_by(species,fleet)%>%
    dplyr::summarise(catch.model = mean(catch.model,na.rm=T))

  catch.mag.ref = catch.ref %>%
    dplyr::group_by(species,fleet)%>%
    dplyr::summarise(catch.ref = mean(ref.value,na.rm=T))


  catch.mag.all = catch.mag.model %>%
    dplyr::left_join(catch.mag.ref)%>%
    dplyr::mutate(catch.rel = catch.ref/catch.model,
                  catch.magnitude = ifelse(catch.rel > (1 - relChangeThreshold) & catch.rel < (1+relChangeThreshold),T,F)
    )

  #distance
  catch.cog.model =catch.fleet %>%
    dplyr::filter(time >= (max.yr - nYrs))%>%
    dplyr::left_join(boxes)%>%
    dplyr::mutate(catch.wgt.x = inside_long * atoutput,
                  catch.wgt.y = inside_lat * atoutput)%>%
    dplyr::group_by(species,fleet,time)%>%
    dplyr::summarise(cog.x = sum(catch.wgt.x,na.rm=T),
                     cog.y = sum(catch.wgt.y,na.rm=T),
                     catch.tot = sum(atoutput,na.rm=T))%>%
    dplyr::mutate(cog.x = cog.x/catch.tot,
           cog.y = cog.y/catch.tot)%>%
    dplyr::group_by(species,fleet)%>%
    dplyr::summarise(cog.x.model = mean(cog.x, na.rm=T),
                     cog.y.model = mean(cog.y, na.rm=T))

  catch.cog.ref =catch.ref %>%
    dplyr::left_join(boxes)%>%
    dplyr::mutate(catch.wgt.x = inside_long * ref.value,
                  catch.wgt.y = inside_lat * ref.value)%>%
    dplyr::group_by(species,fleet)%>%
    dplyr::summarise(cog.x = sum(catch.wgt.x,na.rm=T),
                     cog.y = sum(catch.wgt.y,na.rm=T),
                     catch.tot = sum(ref.value,na.rm=T))%>%
    dplyr::mutate(cog.x.ref = cog.x/catch.tot,
           cog.y.ref = cog.y/catch.tot)

  catch.cog.all = catch.cog.model %>%
    dplyr::left_join(catch.cog.ref)%>%
    dplyr::mutate(dist = sqrt((cog.x.ref - cog.x.model)^2 + (cog.y.ref-cog.y.model)^2),
                  catch.dist = ifelse(dist <= min.dist,T,F))

  #combine for output
  catch.diag.all = catch.mag.all %>%
    dplyr::left_join(catch.cog.all)%>%
    dplyr::select(species,fleet,catch.model,catch.ref, cog.x.model,cog.x.model,cog.x.ref,cog.y.model,cog.y.ref,dist,catch.magnitude,catch.dist)

  return(catch.diag.all)
}
