

#' SecondaryFit class
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' The `SecondaryFit` class contains a secondary inactivation model fitted to a set of
#' parameters obtained under static environmental conditions.
#' Its constructor is [fit_inactivation_secondary()].
#'
#' It is a subclass of list with the items:
#'
#' - ...:...
#'
#' @name SecondaryFit
#'
NULL

#' @describeIn SecondaryFit vector of fitted model parameters.
#'
#' @param object an instance of [SecondaryFit].
#' @param ... ignored
#'
#' @importFrom stats coef
#'
#' @export
#'
coef.SecondaryFit <- function(object, ...) {

  if (object$algorithm == "regression") {

    coef(object$fit_results)

  } else {
    object$fit_results$bestpar
  }

}

#' @describeIn SecondaryFit statistical summary of the fit.
#'
#' @param object Instance of [SecondaryFit]
#' @param ... ignored
#'
#' @export
#'
summary.SecondaryFit <- function(object, ...) {

  out <- summary(object$fit_results)

  # if (object$algorithm != "MCMC") {  # The summary of MCMC is a data.frame, so this would add a column
  #   out$logbase_mu <- object$logbase_mu
  #   out$logbase_logN <- object$logbase_logN
  # }

  out

}
