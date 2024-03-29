% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inactivation_constant.R
\name{predict_constant_inactivation}
\alias{predict_constant_inactivation}
\title{Predicting inactivation under constant environmental conditions}
\usage{
predict_constant_inactivation(times, model_name, model_pars, check = TRUE)
}
\arguments{
\item{times}{a numeric vector of time points where to calculate the prediction}

\item{model_name}{a model identifier as per \code{\link[=primary_model_data]{primary_model_data()}}}

\item{model_pars}{a numeric vector (or list) describing the values of the model
parameters as defined by \code{\link[=primary_model_data]{primary_model_data()}}}

\item{check}{whether to due some model checks. \code{TRUE} by default}
}
\value{
a tibble with two columnes ('time' and 'logN') with the results of the
calculation.
}
\description{
Predictions using only primary inactivation models. In principle, these models
would be valid for constant environmental conditions.
}
