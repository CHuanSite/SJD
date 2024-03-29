% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weightData.R
\name{weightData}
\alias{weightData}
\title{Weighting Data Set}
\usage{
weightData(dataset, weighting)
}
\arguments{
\item{dataset}{A list of data sets}

\item{weighting}{A vector of weighting constant for each data set}
}
\value{
A list of weighted data sets
}
\description{
To weight each data set based on input weighting vector
}
\examples{
dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
weighting = c(1, 2)
weighted_dataset = weightData(dataset, weighting)

}
\keyword{weighting}
