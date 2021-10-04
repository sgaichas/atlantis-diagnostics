#' Test for max fish size
#'
#' Calculate average size of fish through time and check to see if any fish exceed the
#' user supplied observed maximum size. Both maximum weight and maximum length are used
#'
#'@param nc A character string. Path to location of main nc file.
#'@param bgm A character string. Path to location of box coordinates bgm file.
#'@param init A character string. Path to location of initial conditions nc file.
#'@param fgs A character string. Path to location of functional groups file.
#'@param prm_run A character string. Path to location of run parameter file file.
#'@param prm_biol A character string. Path to location of biology parameter file.
#'@param speciesStats Data frame. Must contain at least 3 columns labeled \code{code} - Atlantis species codes
#',\code{maxObsWeight} - a value (g) indicating the maximum weight the species should weigh,
#'\code{maxObsLength} - a value (cm) indicating the maximum length the species should grow to.
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for large fish.
#'(Default = NULL, uses all species)
#'
#'@return Returns a data frame indicating which species meet defined max weight/length criterion.
#'
#'\item{species}{The common name of the species/functional group}
#'\item{code}{Atlantis Code for species/functional group}
#'\item{polygon}{Polygon in which the largest individual occupies}
#'\item{agecl}{Age class of the largest individual}
#'\item{time}{Time in which the largest individual was found (Decimal years)}
#'\item{maxLength}{Maximum length of species in Atlantis (cm)}
#'\item{maxObsLength}{Maximum observed length of species (cm) - from litereature }
#'\item{maxMeanWeight}{Maximum mean weight of indivdual in Atlantis (g)}
#'\item{maxObsWeight}{Maximum observed weight of species (g) - from litereature}
#'\item{passW}{Boolean indicating whether species maximum weight (from Atlantis) falls below max observed weight}
#'\item{ratioW}{ratio of atlantis max weight to observed max weight}
#'\item{passL}{Boolean indicating whether species maximum length (from Atlantis) falls below max observed length}
#'\item{ratioL}{ratio of atlantis max length to observed max length}
#'
#'
#'@importFrom magrittr %>%
#'@importFrom rlang .data
#'
#'@export



diag_maxsize <- function(nc,bgm,init,fgs,prm_run,prm_biol,speciesStats,speciesCodes = NULL){

  # select only the columns required
  speciesStats <- speciesStats %>%
    dplyr::select(.data$code,.data$maxObsLength,.data$maxObsWeight)

  #Get boundary box
  bboxes <-   atlantistools::get_boundary(boxinfo = atlantistools::load_box(bgm))
  #Get epibenthic biopool groups
  bio.pools <-  atlantistools::load_bps(fgs,init)
  #Read in box properties
  vol.dz <-  atlantistools::load_nc_physics(nc = nc,
                                          select_physics = c('volume','dz'),
                                          prm_run = prm_run,
                                          bboxes = bboxes)
  #Get biomass conversion scalar
  bio.conv <-  atlantistools::get_conv_mgnbiot(prm_biol)


  speciesLookup <- atlantistools::load_fgs(fgs) %>%
    dplyr::select(.data$LongName, .data$Code)

  # get group names
  group.names <-  atlantistools::get_groups(fgs)
  # get age groups - 10 cohorts
  groups.age <-  atlantistools::get_age_groups(fgs)
  # get age groups acronyms  - 10 cohorts
  groups.age.acronym <-  atlantistools::get_age_acronyms(fgs)
  # get groups not with 10 cohorts
  groups.bp <-  group.names[!group.names %in% groups.age]

  # list of variables to pull from main nc file.
  # Needed for biomass calculation. Each variable resides in list element
  vars = list('Nums','StructN','ResN','N')
  group.types = list(groups.age,groups.age,groups.age,groups.bp)
  rawdata.main = Map(atlantistools::load_nc,
                     select_variable = vars,
                     select_groups = group.types,
                     MoreArgs = list(nc = nc,
                                     bps = bio.pools,
                                     fgs = fgs,
                                     prm_run = prm_run,
                                     bboxes = bboxes ))

  # calculate biomass for species,age, polygon, layer, time
  spatial.biomass = atlantistools::calculate_biomass_spatial(nums = rawdata.main[[1]],
                                                             sn = rawdata.main[[2]],
                                                             rn = rawdata.main[[3]],
                                                             n = rawdata.main[[4]],
                                                             vol_dz = vol.dz,
                                                             bio_conv = bio.conv,
                                                             bps = bio.pools)
  # grab numbers in time and space
  spatialNumbers = rawdata.main[[1]] %>%
    dplyr::rename(numbers = .data$atoutput)
  # filter biomass for species with 10 cohorts and convert to kilograms
  spatialBiomass <- spatial.biomass %>%
    dplyr::filter(.data$species %in% unique(spatialNumbers$species)) %>%
    dplyr::rename(biomass = .data$atoutput) %>%
    dplyr::mutate(biomass = 1E6*.data$biomass)

  # join numbers with biomass and calculate mean weight of an individivual in age, polygon, layer, time
  jj <- spatialNumbers %>% dplyr::left_join(.,spatialBiomass,by = c("species", "agecl", "polygon", "layer", "time")) %>%
    dplyr::filter(!is.na(.data$biomass)) %>%
    dplyr::mutate(meanWeight = .data$biomass/.data$numbers)

  # select time and space where max occurs
  boxlayer <- jj %>% dplyr::group_by(.data$species,.data$polygon,.data$agecl, .data$time) %>%
    dplyr::summarise(maxMeanWeight = max(.data$meanWeight),.groups="drop")
  maxVal <- boxlayer %>%
    dplyr::group_by(.data$species) %>%
    dplyr::summarise(maxVal = max(.data$maxMeanWeight),.groups="drop")
  boxlayer <- boxlayer %>%
    dplyr::left_join(.,maxVal,by = "species") %>%
    dplyr::filter(.data$maxMeanWeight == .data$maxVal) %>%
    dplyr::select(-.data$maxVal)


  # get length weight params from model W = aL^b
  # note: estimating L from W using this relationship is incorrect
  lenWeightParams <- atlantistools::prm_to_df(prm_biol,fgs,group = groups.age.acronym,parameter = c("li_a","li_b"))

  # find the maximum weight for each species compare to large fish values from speciesStats
  # do the same for the maximum length using Weight-length relationship
  largeFish <- boxlayer %>%
    dplyr::left_join(.,speciesLookup,by=c("species"="LongName")) %>%
    dplyr::rename(code = .data$Code) %>%
    dplyr::left_join(.,speciesStats,by = "code") %>%
    dplyr::mutate(passW = ifelse(.data$maxMeanWeight < .data$maxObsWeight,TRUE,FALSE)) %>%
    #dplyr::select(-.data$scientificName,-.data$Common_Name) %>%
    dplyr::relocate(.data$code,.after = .data$species) %>%
    dplyr::relocate(.data$maxMeanWeight, .before = .data$maxObsWeight) %>%
    dplyr::mutate(ratioW = .data$maxMeanWeight/.data$maxObsWeight) %>%
    dplyr::left_join(.,lenWeightParams,by = "species") %>%
    dplyr::mutate(maxLength = (.data$maxMeanWeight/.data$li_a)^(1/.data$li_b)) %>%
    dplyr::select(-.data$li_a,-.data$li_b) %>%
    dplyr::relocate(.data$maxLength,.after = .data$time) %>%
    dplyr::mutate(passL = ifelse(.data$maxLength < .data$maxObsLength,TRUE,FALSE)) %>%
    dplyr::mutate(ratioL = .data$maxLength/.data$maxObsLength) %>%
    dplyr::arrange(.data$passW,.data$passL)


  # filter codes supplies by user
  if(!is.null(speciesCodes)) {
    largeFish <- largeFish %>%
      dplyr::filter(.data$code %in% speciesCodes)
  }



  return(largeFish)

}
