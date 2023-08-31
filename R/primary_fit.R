
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
#' @param fit_data aa
#' @param model_name description
#' @param start description
#' @param known description
#' @param upper description
#' @param lower description
#' @param check description
#'
#' @importFrom FME modFit modCost
#' 
#' @returns An instance of modFit with the fitted model.
#'
fit_primary <- function(fit_data,
                        model_name,
                        start,
                        known,
                        upper = NULL,
                        lower = NULL,
                        # approach_logN0 = "unique",  # or "logS" or "different",
                        # secondary_models = NULL
                        # algorithm = "regression",
                        # env_conditions = NULL,
                        # niter = NULL,
                        # ...,
                        check = TRUE
) {

  ## Check the model parameters

  if (isTRUE(check)) {

    check_primary_pars(model_name, c(start, known))

  }

  ## Fit the model
  
  if ( !is.null(upper) ) {
    upper <- upper
  } else {
    upper <- Inf
  }
  
  if ( !is.null(lower) ) {
    lower <- lower
  } else {
    lower <- -Inf
  }

  my_fit <- modFit(primary_residuals,
                   unlist(start),
                   fit_data = fit_data,
                   model_name = model_name,
                   known = unlist(known),
                   upper = upper,
                   lower = lower
                   # logbase_logN = logbase_logN,
                   # ...
                   )

  ## Output

  my_fit

}


