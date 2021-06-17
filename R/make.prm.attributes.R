#' Internal
#'
#' Code taken and adapted from shinyrAtlantis
#'
#' @param prmDir Character vector. Path to biology parameter file
#' @param fgs Data frame. Contents of the functional group file. \code{atlantisom::load_fgs}
#'
#' @noRd

make.prm.attributes <- function(prmDir, fgs){
  prm <- readLines(prmDir) # read in the biological parameter file
  grp.att <- fgs %>%
    dplyr::select(Code,Name,GroupType)

  def.grp.file <- system.file("extdata", "grpTemplates.csv", package = "shinyrAtlantis")
  df.prms.all <- read.csv(file = def.grp.file, header = TRUE)
  tmplts <- df.prms.all$Template # group templates to search for
  Codes <- grp.att$Code # groups to search for

  for(tmplt in tmplts) { # look for each template
    cat("-")
    p.vals <- rep(NA, length(Codes))
    i <- 0
    for (xxx in Codes) { # look for each Code
      i <- i + 1 # xxx index
      txt.find <- gsub(pattern = "XXX", replacement = xxx, x = tmplt)
      j <- grep(pattern = txt.find, x = prm, value = FALSE) # file row(s)
      if (length(j) > 0) { # found the parameter
        jnew <- NULL
        for (jj in 1:length(j)) {
          # Valid row is when tmplt is the first entry and second is a number
          text.split <- unlist(stringr::str_split(
            gsub(pattern = "[ \t]+", x = prm[j[jj]], replacement = " "), " "))
          if (text.split[1] == txt.find) {
            jnew <- c(jnew,j[jj]) # add the row that satisfies the criteria
          }
        }
        j <- jnew # use this list of rows as they are valid
        if (length(j) == 1) { # a single row is found
          # get parameter value (1 after nums)
          p.val <- as.numeric(unlist(stringr::str_split(prm[j],"[\t ]+"))[2])
          p.vals[i] <- p.val
        }
      }
    }
    grp.att$tmplt <- p.vals # add the group attribute column
    names(grp.att)[names(grp.att) == "tmplt"] <- as.character(tmplt) # rename column
  }
  cat("\n")

  return(grp.att)
}
