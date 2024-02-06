

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

  out

}

#' @describeIn SecondaryFit vector of model residuals.
#'
#' @param object Instance of [SecondaryFit]
#' @param ... ignored
#'
#' @importFrom stats residuals
#'
#' @export
#'
residuals.SecondaryFit <- function(object, ...) {
  
  if (object$algorithm == "MCMC") {
    
    pred <- predict(object)
    
    pred - object$data$logN
    
  } else {
    
    residuals(object$fit_results)
    
  }
  
}

#' @describeIn SecondaryFit variance-covariance matrix of the model, estimated
#' as 1/(0.5*Hessian) for regression and as the variance-covariance of the draws
#' for MCMC
#'
#' @param object an instance of [SecondaryFit]
#' @param ... ignored
#' 
#' @importFrom stats cov
#'
#' @export
#'
vcov.SecondaryFit <- function(object, ...) {
  
  if (object$algorithm == "MCMC") {
    
    cov(object$fit_results$pars)
    
  } else {
    
    # The code has been adapted from the one of summary.modFit
    
    covar  <- try(solve(0.5*object$fit_results$hessian), silent = TRUE)
    
    if (!is.numeric(covar)) {
      warning("Cannot estimate covariance; system is singular")
      
      param  <- object$par
      p      <- length(param)
      
      covar <- matrix(data = NA, nrow = p, ncol = p)
    }
    
    covar
    
  }
}

#' @describeIn SecondaryFit deviance of the model.
#'
#' @param object an instance of [SecondaryFit]
#' @param ... ignored
#'
#' @importFrom stats deviance
#'
#' @export
#'
deviance.SecondaryFit <- function(object, ...) {
  
  if (object$algorithm == "MCMC") {
    
    sum(residuals(object)^2)
    
  } else {
    deviance(object$fit_results)
  }
  
}

#' @describeIn SecondaryFit loglikelihood of the model
#'
#' @param object an instance of SecondaryFit
#' @param ... ignored
#'
#' @export
#'
logLik.SecondaryFit <- function(object, ...) {
  
  if (object$algorithm == "regression") {
    
    n <- nrow(object$data)
    sigma <- sqrt(object$fit_results$ssr/object$fit_results$df.residual)
    
    lL <- - n/2*log(2*pi) -n/2 * log(sigma^2) - 1/2/sigma^2*object$fit_results$ssr
    
    lL
    
  } else {
    
    n <- nrow(object$data)
    SS <- min(object$fit_results$SS, na.rm = TRUE)
    
    df <- n - length(coef(object))
    
    sigma <- sqrt(SS/df)
    
    lL <- - n/2*log(2*pi) -n/2 * log(sigma^2) - 1/2/sigma^2*SS
    
    lL
    
  }
}

#' @describeIn SecondaryFit Akaike Information Criterion
#'
#' @param object an instance of SecondaryFit
#' @param ... ignored
#' @param k penalty for the parameters (k=2 by default)
#'
#' @importFrom stats logLik
#'
#' @export
#'
AIC.SecondaryFit <- function(object, ..., k=2) {
  
  ## Normal AIC
  
  p <- length(coef(object))
  
  lL <- logLik(object)
  
  AIC <- 2*p - 2*lL
  
  ## Calculate the penalty
  
  n <- nrow(object$data)
  
  penalty <- (k*p^2 + k*p)/(n - p - 1)
  
  ## Return
  
  AIC + penalty
  
}

#' @describeIn SecondaryFit vector of fitted values.
#'
#' @param object an instance of [SecondaryFit]
#' @param ... ignored
#'
#' @export
#'
fitted.SecondaryFit <- function(object, ...) {
  
  res <- secondary_residuals(coef(object), 
                             object$model_name, 
                             object$data, 
                             known = object$known) 
  
  object$data$my_par - res
  
}


#' @describeIn SecondaryFit vector of model predictions.
#' 
#' @importFrom stats fitted
#'
#' @param object an instance of [SecondaryFit]
#' @param ... ignored
#' @param newdata tibble describing the environmental conditions as in [fit_inactivation_secondary()].
#' If `NULL` (default), uses the environmental condition of the fitting.
#'
#' @export
#'
predict.SecondaryFit <- function(object, newdata = NULL, ...) {
  
  if ( is.null(newdata) ) {
    fitted(object)
  } else {
    
    fake <- newdata %>% mutate(my_par = 1)
    
    res <- secondary_residuals(coef(object), 
                               object$model_name, 
                               fake, 
                               known = object$known) 
    
    fake$my_par - res
    
  }
  
}

#' @describeIn SecondaryFit compares the fitted model against the data.
#' 
#' @param x an instance of SecondaryFit
#' @param y ignored
#' @param ... ignored
#' @param type type of plot to make. Either 1 (observed vs predicted) or 2 (predictived vs residuals).
#'
#' @importFrom stats fitted
#' @export
#'
plot.SecondaryFit <- function(x, y=NULL, ..., type = 1) {

  d <- x$data %>%
    mutate(res = residuals(x),
           pred = fitted(x)
           )
  
  if (type == 1) {
    
    ggplot(d, aes(x = .data$my_par, y = .data$pred)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      geom_abline(slope = 1, intercept = 0, linetype = 2)
    
  } else if (type == 2) {
    
    ggplot(d, aes(x = .data$pred, y = .data$res)) +
      geom_point() +
      geom_smooth(se = FALSE) +
      geom_hline(yintercept = 0, linetype = 2) 
    
  } else {
    stop(paste("Type must be 1 or 2, got:", type))
  }
  
}









