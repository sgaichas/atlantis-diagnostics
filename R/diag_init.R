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
    dplyr::filter(IsTurnedOn == 1) %>%
    dplyr::select(Name) %>%
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
