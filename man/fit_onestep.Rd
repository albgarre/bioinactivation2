% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/one_step_fit.R
\name{fit_onestep}
\alias{fit_onestep}
\title{One-step fitting of isothermal inactivation data}
\usage{
fit_onestep(
  fit_data,
  model_name,
  start,
  known,
  upper = NULL,
  lower = NULL,
  secondary_models,
  algorithm,
  niter = NULL,
  ...
)
}
\arguments{
\item{fit_data}{a tibble (or data.frame) with the data for the fitting.}

\item{model_name}{a model identifier according to \code{\link[=primary_model_data]{primary_model_data()}}}

\item{start}{a named numeric vector of initial guesses for the parameters as
defined in \code{\link[=fit_inactivation]{fit_inactivation()}}}

\item{known}{a named numeric vector of known parameters}

\item{upper}{a named numeric vectors of upper bounds for the parameters}

\item{lower}{a named numeric vector of lower bounds for the parameters}

\item{secondary_models}{a nested list defining the secondary models for each
parameter as defined in \code{\link[=fit_inactivation]{fit_inactivation()}}}

\item{algorithm}{one of \code{"regression"} (default) or \code{"MCMC"}}

\item{niter}{number of MC iterations. Ignored if \verb{algorithm == "regression}}

\item{...}{additional arguments for \code{\link[=modFit]{modFit()}} or \code{\link[=modMCMC]{modMCMC()}}}
}
\value{
An instance of \link{InactivationFit}
}
\description{
One-step fitting of isothermal inactivation data
}
