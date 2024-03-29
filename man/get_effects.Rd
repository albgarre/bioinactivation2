% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_effects}
\alias{get_effects}
\title{Effects of the environmental conditions on the reference value}
\usage{
get_effects(this_par, env_conditions)
}
\arguments{
\item{this_par}{a nested list defining the secondary model as in \code{\link[=fit_inactivation]{fit_inactivation()}}.
Note that anti-log transformations must be applied before calling this function.}

\item{env_conditions}{a named numeric vector with the values fo the environmental conditions
(just 1 time point)}
}
\value{
a named numeric vector with the effect of each factor on the reference
value.
}
\description{
The function calculates the effect of each environmental factor on the reference
value according to the secondary model defined.
}
