demo.part <- function(partno){
  guides <- c("We want to demonstrate a little and easy demo session.\nFirst of all, we have to include the library and to load the provided\nexpression data on acute lymphocytic leukemia (ALL).",
              "We now have an data object called 'stjude'.\nLet's have a look at the data.",
              "Among the components of this list are the identifiers of the genes and the labels\nof the samples in the dataset.",
              "There are ten different classes of tumor patients.\nThe next questions is: how many probe sets are on chromosome 1?",
              "OK, now lets plot the expression values on the sixth chromosome of the third sample.\nFor smoothing, we use a 'radial basis function (rbf)' kernel:",
              "Now, we want to look for significantly differentially expressed chromosome regions.\nWe investigate chromosome 6 for interesting regions with patient class 'T',\nwhich means 'T-lymphocyte ALL'. We test the class 'T' against the common other\nclasses defined as background in this framework. First we try to find the\noptimal parameter settings for the smoothing kernel function, for which we chose\na k-nearest-neighbor kernel.",
              "Have a look at the average residual sum-of-squares for different\nsettings for parameter k = number of neighbors.",
              "Now, we use the 'evalScoring' function to search for interesting regions on\nchromosome 6.\nWe use the k-nearest-neighbor kernel with the optimal parameter settings,\nestimated above. Since we did cross-validation above, we can omit it now.\nThe process may take some time due to the number of permutations.",
              "Have a look at the result.",
              "We can also use the plot function to generate an HTML page. For this, we\nhave to set the argument 'output' to 'html'.See the result in your browser.\nWe get some information about genes in highlighted (=significant) regions.")
  cat(guides[partno],fill=TRUE)
  cat("<RETURN> to continue...")
  invisible(readline())
}#demo.part
