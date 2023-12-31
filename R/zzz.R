.onLoad <- function(libname, pkgname) {
  ## nothing to do here for the moment
}

.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.19")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
    
    packageStartupMessage("Loading MicroArray Chromosome Analysis Tool...\n",
                        "Loading required packages...\n")
  
  if (!("stjudem" %in% .packages(all.available=TRUE))){
    packageStartupMessage(
      "You need package 'stjudem' if you want to see the demo and examples.\n",
      "Type 'loaddatapkg(\"stjudem\")' for automatic install!\n\n")
  } # if data package not installed
  packageStartupMessage("Type 'demo(macatdemo)' for a quick tour...\n")
  ## show vignette in windows menu  
  if(.Platform$OS.type=="windows" && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("macat")
  }
} # .onAttach

.onUnload <- function( libpath ) {
  ## nothing to do here for the moment
} # .onUnload
