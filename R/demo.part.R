demo.part <- function(partno){
  guides <- c("We want to demonstrate a little and easy demo session.\nFirst of all, we have to include the library and to load the provided expression data\non acute lymphocytic leukemia (ALL).",
              "We now have an data object called 'stjude'.\nLet's have a look at the data.",
              "Among the components of this list are the identifiers of the genes and the labels of the samples in the dataset.",
              "There are ten different classes of tumor patients.\nThe next questions is: how many probe sets are on chromosome 1?",
              "OK, now lets plot the expression values on the sixth chromosome of the third sample.\nFor smoothing, we use a 'radial basis function (rbf)' kernel:",
              "Now we want to look for chromosomal aberrations.\nWe look again on chromosome 6 for some specific aberrations for patient class 'T',\nwhich means 'T-lymphocyte ALL'.\nWe test the class 'T' against the common other classes defined as background in this framework.\nFirst we use the 'evalScoring' function to build an MACATevalScoring object.\nThis may take some time due to the number of permutations.\nThis time we use the default rbf kernel to smooth the scores.",
              "Now let's use a k-nearest-neighbor kernel:",
              "We can compare the two results by using the plot function:",
              "Now we can use the plot function to generate an HTML-page on-the-fly.\nTherefore we have to set the parameter output to 'html'.\nSee the result in your browser.\nWe get some useful information about genes in highlighted (=significant) regions.")
  cat(guides[partno],fill=TRUE)
  cat("<RETURN> to continue...")
  invisible(readline())
}#demo.part
