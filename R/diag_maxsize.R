#' Test for max fish size
#'
#' Calculate average size of fish through time and check to see if any fish exceed the
#' user supplied observed maximum size.
#'
#'@param nc A character string. Path to location of model run output files.
#'@param bgm A character string. Prefix added to all atlantis output files.
#'@param init A character string. Prefix added to all atlantis output files.
#'@param fgs A character string. Prefix added to all atlantis output files.
#'@param prm_run A character string. Prefix added to all atlantis output files.
#'@param prm_biol A character string. Prefix added to all atlantis output files.#'
#'@param speciesStats Data frame.
#'@param speciesCodes Character vector. A vector of Atlantis species codes in which to test for persistence.
#'(Default = NULL, uses all species)
#'
#'@importFrom magrittr %>%
#'@importFrom rlang .data
#'
#'
#'@return Returns a data frame of species which do not meet defined max size criterion.
#'
#'\item{species}{The common name of the species/functional group}
#'\item{code}{Atlantis Code for species/functional group}
#'
#'@export
#'
#'
#'@examples
#'\dontrun{
#'
#'}


diag_maxsize <- function(nc,bgm,init,fgs,prm_run,prm_biol,speciesStats,speciesCodes = NULL){

  #Get boundary box
  bboxes =  atlantistools::get_boundary(boxinfo = atlantistools::load_box(bgm))
  #Get epibenthic biopool groups
  bio.pools = atlantistools::load_bps(fgs,init)
  #Read in box properties
  vol.dz = atlantistools::load_nc_physics(nc = nc,
                                          select_physics = c('volume','dz'),
                                          prm_run = prm_run,
                                          bboxes = bboxes)
  #Get biomass conversion scalar
  bio.conv = atlantistools::get_conv_mgnbiot(prm_biol)


  speciesLookup <- atlantistools::load_fgs(fgs) %>%
    dplyr::select(.data$LongName, .data$Code) %>%
    tibble::as_tibble()

  # get group names
  group.names = atlantistools::get_groups(fgs)
  # get age groups - 10 cohorts
  groups.age = atlantistools::get_age_groups(fgs)
  # get groups not with 10 cohorts
  groups.bp = group.names[!group.names %in% groups.age]

  # list of variables ot pull from main nc file.
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
    dplyr::rename(numbers = atoutput)
  # filter biomass for species with 10 cohorts and convert to kilograms
  spatialBiomass <- spatial.biomass %>%
    dplyr::filter(.data$species %in% unique(spatialNumbers$species)) %>%
    dplyr::rename(biomass = .data$atoutput) %>%
    dplyr::mutate(biomass = 1E6*.data$biomass)

  # join numbers with biomass and calculate mean weight of an individivual in age, polygon, layer, time
  jj <- spatialNumbers %>% dplyr::left_join(.,spatialBiomass,by = c("species", "agecl", "polygon", "layer", "time")) %>%
    dplyr::filter(!is.na(.data$biomass)) %>%
    dplyr::mutate(meanWeight = .data$biomass/.data$numbers)

  # select time and spac where max occurs
  boxlayer <- jj %>% dplyr::group_by(.data$species,.data$polygon, .data$time) %>%
    dplyr::summarise(maxMeanWeight = max(.data$meanWeight),.groups="drop")
  maxVal <- boxlayer %>%
    dplyr::group_by(.data$species) %>%
    dplyr::summarise(maxVal = max(.data$maxMeanWeight),.groups="drop")
  boxlayer <- boxlayer %>%
    dplyr::left_join(.,maxVal,by = "species") %>%
    dplyr::filter(maxMeanWeight == .data$maxVal) %>%
    dplyr::select(-.data$maxVal)


  # find the maximum weight for each species compare to large fish values from speciesStats
  largeFish <- boxlayer %>%
    dplyr::left_join(.,speciesLookup,by=c("species"="LongName")) %>%
    dplyr::left_join(.,speciesStats,by = c("Code"="code")) %>%
    dplyr::mutate(pass = ifelse(.data$maxMeanWeight<.data$maxObsWeight,TRUE,FALSE)) %>%
    dplyr::arrange(.data$pass)  %>%
    dplyr::select(-.data$scientificName,-.data$Common_Name) %>%
    dplyr::relocate(.data$Code,.after = .data$species) %>%
    dplyr::relocate(.data$maxMeanWeight, .before = .data$maxObsWeight) %>%
    dplyr::mutate(ratio = .data$maxMeanWeight/.data$maxObsWeight)

  # filter codes supplies by user
  if(!is.null(speciesCodes)) {
    largeFish <- largeFish %>%
      dplyr::filter(Code %in% speciesCodes)
  }



  return(largeFish)

}
