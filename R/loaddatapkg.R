loaddatapkg <- function(mydatapkg,installDir=.libPaths()[1]){
  if (mydatapkg %in% .packages(all.available=TRUE)) {
    library(mydatapkg,character.only=TRUE)
  } else {
    if ("reposTools" %in% .packages(all.available=TRUE)) {
      library(reposTools)
      myrepos <- getReposEntry("http://compdiag.molgen.mpg.de/software/Rrepos")
      insttry <- try(install.packages2(mydatapkg,myrepos,lib=installDir,type="Source",develOK=TRUE),silent=TRUE)      
      libtry <- try(library(mydatapkg,character.only=TRUE),silent=TRUE)
      if ((class(insttry)=="try-error")|(class(libtry)=="try-error")){
        stop(paste("Couldn't load data package '",mydatapkg,"'!\nPlease, try to obtain it from\n\nhttp://compdiag.molgen.mpg.de/software/macat.shtml\n\nand install it yourself.\n",sep=""),call.=FALSE)
      } # if (class(insttry)=="try-error")
    } else { # no reposTools
      stop(paste("Need Package 'reposTools' for automatic download. !\nPlease, try to obtain package '",mydatapkg,"' manually from \n\nhttp://compdiag.molgen.mpg.de/software\n\n install it and rerun this demo.\n",sep=""),call.=FALSE)
    } # else if no reposTools
  } # else package not available
} # loaddatapkg


