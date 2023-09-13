
#' Residuals for fitting primary inactivation models
#' 
#' The function is prepared to be called by modFit or modMCMC
#' 
#' @param this_p named numeric vector with the candidate model parameters
#' @param fit_data data for the fit, as a tibble (or data frame) with two columns: time and logN
#' @param model_name model identifier
#' @param known named numeric vector with the known parameters
#' 
#' @returns an instance of modCost
#'
primary_residuals <- function(this_p,
                              fit_data,
                              model_name,
                              known
                              # logbase_logN = 10
) {

  ## Make the prediction

  times <- sort(unique(fit_data$time))

  pars <- c(this_p, known)
  my_model <- as.list(pars)
  my_model$model <- model_name

  pred <- predict_inactivation(times,
                       my_model,
                       check = FALSE
                       # logbase_logN = logbase_logN
                       )

  ## Return the residuals

  modCost(model = as.data.frame(pred$simulation),
          obs = as.data.frame(fit_data)
          )

}


#' Fitting of primary inactivation models
#' 
#' @param fit_data tibble (or data.frame) with the data used for fitting. It must have 
#' a column named `time` and one named `logN`
#' @param model_name a model identifier according to [primary_model_data()]
#' @param start a named numeric vector of initial guesses for the model parameters. They
#' must be named according to [primary_model_data()]
#' @param known named numeric vector of known parmeters
#' @param upper named numeric vector of upper limits for the parameters
#' @param lower named numeric vector of lower limits for the parameters
#' @param check whether to do some basic checks of the model parameters
#' @param ... additional arguments passed to [modFME()] or [modFit()]
#' @param algorithm fitting algorithm: "regression" (default) or "MCMC"
#' @param niter number of MCMC iterations. Ignored if `algorithm == "regression"`
#'
#' @importFrom FME modFit modCost
#' 
#' @returns An instance of [InactivationFit] with the fitted model.
#'
fit_primary <- function(fit_data,
                        model_name,
                        start,
                        known,
                        upper = NULL,
                        lower = NULL,
                        algorithm = "regression",
                        niter = NULL,
                        ...,
                        check = TRUE
) {

  ## Check the model parameters

  if (isTRUE(check)) {

    check_primary_pars(model_name, c(start, known))

  }

  ## Set up the bounds
  
  if (is.null(upper)) upper <- Inf
  if (is.null(lower)) lower <- -Inf
  
  ## Fit the model
  
  if (algorithm == "regression") {
    
    my_fit <- modFit(primary_residuals,
                     unlist(start),
                     fit_data = fit_data,
                     model_name = model_name,
                     known = unlist(known),
                     upper = upper,
                     lower = lower,
                     ...
                     )
    
    ## Calculate the best prediction
    
    p <- c(coef(my_fit), unlist(known))
    primary_model <- as.list(p)
    primary_model$model <- model_name
    t <- seq(0, max(fit_data$time), length = 1000)
    
    best_prediction <- predict_inactivation(t, primary_model)
    
    ## Prepare the output
    
    out <- list(
      approach = "primary",
      algorithm = "regression",
      data = fit_data,
      guess = start,
      known = known,
      primary_model = model_name,
      fit_results = my_fit,
      best_prediction = best_prediction,
      sec_models = NULL,
      env_conditions = NULL,
      niter = NULL
    )
    
    class(out) <- c("InactivationFit", class(out))
    
    
  } else if (algorithm == "MCMC") {
    
    my_fit <- modMCMC(primary_residuals,
                     unlist(start),
                     fit_data = fit_data,
                     model_name = model_name,
                     known = unlist(known),
                     upper = upper,
                     lower = lower,
                     niter = niter,
                     ...
    )
    
    ## Calculate the best prediction
    
    p <- c(my_fit$bestpar, unlist(known))
    primary_model <- as.list(p)
    primary_model$model <- model_name
    t <- seq(0, max(fit_data$time), length = 1000)
    
    best_prediction <- predict_inactivation(t, primary_model)
    
    ## Prepare the output
    
    out <- list(
      approach = "primary",
      algorithm = "MCMC",
      data = fit_data,
      guess = start,
      known = known,
      primary_model = model_name,
      fit_results = my_fit,
      best_prediction = best_prediction,
      sec_models = NULL,
      env_conditions = NULL,
      niter = NULL
    )
    
    class(out) <- c("InactivationFit", class(out))
    
    
  } else {
    stop(paste("algorithm must be either 'regression' or 'MCMC', got:", algorithm))
  }

  ## Output

  out

}


