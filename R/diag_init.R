#'Initialize Atlantis diagnostics
#'
#'@description \code{diag_init} pulls necessary Atlantis files from google drive and loads
#'them into R objects for input into each diagnostic function.
#'
#'@param config_file file specifying google directory with atlantis model outputs and a local
#'directory for download, and identifying key input files. See example config_file below.
#'
#'@return Returns a diaglist list object containing dataframes and lists:
#' \itemize{
#'  \item{fgs, dataframe of species group characteristics}
#'  \item{fgs.names, vector of species names}
#'  \item{atBtxt, dataframe output of \code{atlantisom::load_bioind}}
#'  \item{atCtxt, dataframe output of \code{atlantisom::load_catch}}
#'  \item{runpar, list of run parameters output of \code{atlantisom::load_runprm}}
#'  \item{noutsteps, integer number of output steps in the model run}
#'  \item{stepperyr, integer number of putput steps per year}
#'  \item{timeall, vector of all model timesteps}
#'  \item{fstepperyr, integer fishery number output steps per year}
#'  \item{boxpars, dataframe of spatial parameters}
#'  \item{boxall, vector of all model box numbers}
#'  \item{YOY, dataframe young of year output}
#'  \item{biol, list of biological parameters output of \code{atlantisom::load_biolprm}}
#' }
#'@export
#'
#'@author Sarah Gaichas
#'
#'@examples
#'\dontrun{
#'# Example config file DiagnosticConfig.R:
#'#------------------------------------
#'# CHECK ALL FILENAMES FOR YOUR MODEL
#'
#'# where are files on the atlantis google drive?
#'g.name <- "Testing/OutForSarah"
#'
#'# which local directory should they write to?
#'d.name <- here::here("diagnostics", "temp")
#'
#'# input file names
#'functional.groups.file <- "neus_groups.csv"
#'biomass.pools.file <- "neus_init.nc"
#'box.file <- "neus_tmerc_RM2.bgm"
#'initial.conditions.file <- "neus_init.nc"
#'fisheries.file <- "neus_fisheries.csv"
#'biol.prm.file <- "at_biology.prm"
#'
#'# define the hindcast period--these are survey years
#'hindcast <- c(1980:2010)
#'#------------------------------------
#'
#'#Usage
#'testdiag <- diag_init("DiagnosticConfig.R")
#'}
#'
#'
diag_init <- function(config_file){

  source(config_file)

  # pull all files from google
  atlantisdrive::pull_from_drive(d.name,
                                 fileList = NULL,
                                 googledriveFolder = g.name)

  # get scenario name from a file inside folder: SSB is unique
  scenario.name <- gsub("SSB.txt","",list.files(path=d.name, pattern = "*SSB.txt"))

  # output file names
  run.prm.file <- list.files(path=d.name, pattern = "*at_run*.xml")
  bioind.file <- paste0(scenario.name, "BiomIndx.txt")
  catch.file <- paste0(scenario.name, "Catch.txt")


  #Load functional groups
  fgs <- atlantisom::load_fgs(dir=d.name,
                              file_fgs = functional.groups.file)
  #Get just the names of active functional groups
  fgs.names <- fgs %>%
    dplyr::filter(.data$IsTurnedOn == 1) %>%
    dplyr::select(.data$Name) %>%
    .$Name

  # should return all model areas
  boxpars <- atlantisom::load_box(d.name, box.file)
  boxall <- c(0:(boxpars$nbox - 1))

  # generalized timesteps all models
  runpar <- atlantisom::load_runprm(d.name, run.prm.file)
  noutsteps <- runpar$tstop/runpar$outputstep
  stepperyr <- if(runpar$outputstepunit=="days") 365/runpar$toutinc
  timeall <- c(0:noutsteps)

  # learned the hard way this can be different from ecosystem outputs
  fstepperyr <- if(runpar$outputstepunit=="days") 365/runpar$toutfinc

  # load the biomass index results
  atBtxt <- atlantisom::load_bioind(d.name, bioind.file, fgs)

  #load the catch results
  atCtxt <- atlantisom::load_catch(d.name, catch.file, fgs)

  # load YOY
  YOY <- atlantisom::load_yoy(d.name, paste0(scenario.name, "YOY.txt"))

  # load biol_prm
  biol <- atlantisom::load_biolprm(d.name, biol.prm.file)


  diaglist <-list("fgs" = fgs,
                  "fgs.names" = fgs.names,
                  "atBtxt" = atBtxt,
                  "atCtxt" = atCtxt,
                  "runpar" = runpar,
                  "noutsteps" = noutsteps,
                  "stepperyr" = stepperyr,
                  "timeall" = timeall,
                  "fstepperyr" = fstepperyr,
                  "boxpars" = boxpars,
                  "boxall" = boxall,
                  "YOY" = YOY,
                  "biol" = biol)

  return(diaglist)

}
