loaddatapkg <- function(mydatapkg,installDir=.libPaths()[1]){
  if (mydatapkg %in% .packages(all.available=TRUE)) {
    library(mydatapkg,character.only=TRUE)
  } else {
    install.packages(mydatapkg, lib=installDir, contriburl=contrib.url("http://www.ebi.ac.uk/huber-srv/data/Rrepos"))
    library(mydatapkg,character.only=TRUE)
  } # else package not available
} # loaddatapkg
