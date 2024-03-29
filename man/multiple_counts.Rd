% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{multiple_counts}
\alias{multiple_counts}
\title{Data for global fitting}
\format{
A list of length two, with each element containing a tibble of 20 rows and 2 columns:
\describe{
\item{time}{treatment time in min}
\item{logN}{observed microbial concentration (log CFU/mL)}
}
}
\usage{
multiple_counts
}
\description{
An example dataset to illustrate \code{\link[=fit_inactivation]{fit_inactivation()}} using the global approach.
It simulates the change in the microbial concentration observed in 2 independent experiments
at dynamic conditions. It should be used in combination with \link{multiple_environments}
}
\keyword{datasets}
