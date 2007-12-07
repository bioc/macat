\name{Python data format}
\alias{pydata}
\alias{pythondata}
\docType{data}
\title{Flat file format}
\description{
This is a description of the text files, usually written by MACAT functions
such as \code{\link{kernelizeToPython}}. These text files can be read into
PYTHON.
Each flat files contains a list of matrices. The outer list has one entry for
each class in the data set, the inner matrices consist of all the sequences from
the specified chromosome belonging to that class. Each row represents one sample,
each column a gene from the respective chromosome.

The lists are delimited by '[' and ']' brackets. The matrices are written as
lists of lists (e.g. [[],[],...,[]] ). 
}

\examples{
# The files look like this:
# [ [[1.11, 1.32, 0.92, ...],[...],[...]],
#   [[0.45, 0.91, 1.84, ...],[...],[...]],
#   [[1.06, 1.59, 0.73, ...],[...],[...]],
# ]  
}
\keyword{internal}

