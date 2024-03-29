% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_inactivation_secondary.R
\name{secondary_residuals}
\alias{secondary_residuals}
\title{Residuals for fitting secondary models}
\usage{
secondary_residuals(this_p, model_name, fit_data, known, output = "vector")
}
\arguments{
\item{this_p}{named numeric vector of parameter values}

\item{model_name}{identifier of the secondary model as per \code{\link[=secondary_model_data]{secondary_model_data()}}}

\item{fit_data}{a tibble (or data.frame) defining the data as in \code{\link[=fit_inactivation_secondary]{fit_inactivation_secondary()}}}

\item{known}{a named numeric vector of known model parameters}

\item{output}{output format. Either 'vector' (default) or 'loglik'}
}
\value{
a numeric vector of model residuals
}
\description{
The function is defined to be called within modFit or modMCMC
}
