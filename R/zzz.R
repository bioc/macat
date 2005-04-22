.First.lib <- function(libname, pkgname){
  cat("Loading MicroArray Chromosome Analysis Tool...\n")
  cat("Loading required packages...\n")
  require(Biobase,quietly=TRUE)
  require(annotate,quietly=TRUE)
  if (!("stjudem" %in% .packages(all.available=TRUE))){
    cat("You need package 'stjudem', if you want to see the demo and examples.\nType 'loaddatapkg(\"stjudem\")' for automatic install!\n\n")
  } # if data package not installed
  cat("Type 'demo(macatdemo)' for a quick tour...\n")
  
  if(.Platform$OS.type=="windows" && require(Biobase) && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("macat")
  }
} # .First.lib

