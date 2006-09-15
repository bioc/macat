loaddatapkg <- function(mydatapkg,installDir=.libPaths()[1]){
  if (mydatapkg %in% .packages(all.available=TRUE)) {
    require(mydatapkg,character.only=TRUE)
  } else {
    install.packages(mydatapkg, lib=installDir, contriburl=contrib.url("http://www.ebi.ac.uk/huber-srv/data/Rrepos"))
    require(mydatapkg,character.only=TRUE)
  } # else package not available
} # loaddatapkg
