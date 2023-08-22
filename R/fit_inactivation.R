


#' Fits an inactivation model using different methods from predictive microbiology
#' 
#' @param method method for model fitting. One of 'primary', 'two-steps', 'one-step',
#' 'dynamic' or 'global'.
#' @param fit_data a tibble with the data to use for model fitting. See sections below for further information.
#' @param primary_model_name aa
#' @param guess description
#' @param known aa
#' @param upper aa
#' @param lower description
#' @param secondary_models description
#' @param algorithm description
#' @param env_conditinos description
#' @param niter description
#' @param check aa
#' 
#' 
#' @export
#'
fit_inactivation <- function(method,
                             fit_data,
                             primary_model_name,
                             guess,
                             known,
                             upper = NULL,
                             lower = NULL,
                             # approach_logN0 = "unique",  # or "logS" or "different",
                             secondary_models = NULL,
                             algorithm = "regression",
                             env_conditions = NULL,
                             niter = NULL,
                             # ...,
                             check = TRUE,
                             # logbase_mu = logbase_logN,
                             # logbase_logN = 10,  # TODO
                             formula = logN ~ time
                             ) {

  if (method == "primary") {
    
    ## Apply the formula
    
    if (length(get.vars(formula)) > 2) {
      stop("Only formulas with 2 terms are supported.")
    }
    
    y_col <- lhs(formula)
    x_col <- rhs(formula)
    
    fit_data <- select(fit_data,
                       time = x_col,
                       logN = y_col
    )

    ## Fit the model

    my_fit <-  fit_primary(fit_data, 
                           primary_model_name, 
                           guess, known,
                           check = check
                           )
    

    ## Calculate the best prediction

    p <- c(coef(my_fit), unlist(known))
    primary_model <- as.list(p)
    primary_model$model <- primary_model_name
    t <- seq(0, max(fit_data$time), length = 1000)

    best_prediction <- predict_inactivation(t, primary_model)

    ## Prepare the output

    out <- list(
      method = method,
      algorithm = "regression",
      data = fit_data,
      guess = guess,
      known = known,
      primary_model = primary_model_name,
      fit_results = my_fit,
      best_prediction = best_prediction,
      sec_models = NULL,
      env_conditions = NULL,
      niter = NULL,
      logbase_logN = NULL
      # approach_logN0 = NULL
    )

    class(out) <- c("InactivationFit", class(out))

    ## Return

    out

  } else if (method == "two-steps") {

    out <- fit_two_step(
      fit_data = fit_data,
      model_name = primary_model_name,
      start = guess,
      known = known,
      upper = upper,
      lower = lower,
      secondary_models = secondary_models
      # approach_logN0 = approach_logN0
    )

    ## Prepare the output

    out$method <- method
    # out$best_prediction <- NA
    # env_conditions <- NULL
    # niter <- NULL
    # logbase_logN = NULL
    # approach_logN0 = NULL

    class(out) <- c("InactivationFit", class(out))

    ## Return

    out



  } else if (method == "one-step") {

    ## Fit the model

    my_fit <- fit_onestep(
      fit_data = fit_data,
      model_name = primary_model_name,
      start = guess,
      known = known,
      upper = upper,
      lower = lower,
      secondary_models = secondary_models,
      algorithm = algorithm,
      niter = niter
      # approach_logN0 = approach_logN0
    )

    ## Extract the secondary and primary models

    p <- c(coef(my_fit), unlist(known))

    primary_model <- list(model = primary_model_name)

    initial <- p[str_detect(names(p), "N0")]
    aa <- initial
    names(aa) <- NULL
    primary_model[[names(initial)]] <- aa

    if (str_detect(primary_model_name, "Geeraerd")) {
      initial <- p[str_detect(names(p), "C0")]
      aa <- initial
      names(aa) <- NULL
      primary_model[[names(initial)]] <- aa
    }

    sec <- convert_dynamic_guess(sec_models, guess, known)

    ## Prepare the output

    out <- list(
      method = method,
      algorithm = "regression",
      data = fit_data,
      guess = guess,
      known = known,
      primary_model = primary_model_name,
      fit_results = my_fit,
      best_prediction = NA,
      sec_models = sec,
      env_conditions = NA,
      niter = niter
      # logbase_logN = NULL
      # approach_logN0 = NULL
    )

    class(out) <- c("InactivationFit", class(out))

    ## Return

    out


  } else if (method == "dynamic") {

    ## Fit the model

    my_fit <- fit_dynamic_inactivation(
      fit_data = fit_data,
      model_name = primary_model_name,
      start = guess,
      known = known,
      upper = upper,
      lower = lower,
      secondary_models = secondary_models,
      algorithm = algorithm,
      env_conditions = env_conditions,
      niter = niter
    )

    ## Calculate the best prediction

    p <- c(coef(my_fit), unlist(known))

    primary_model <- list(model = primary_model_name)

    initial <- p[str_detect(names(p), "N0")]
    aa <- initial
    names(aa) <- NULL
    primary_model[[names(initial)]] <- aa

    if (str_detect(primary_model_name, "Geeraerd")) {
      initial <- p[str_detect(names(p), "C0")]
      aa <- initial
      names(aa) <- NULL
      primary_model[[names(initial)]] <- aa
    }

    t <- seq(0, max(fit_data$time), length = 1000)


    sec <- convert_dynamic_guess(sec_models, guess, known)

    best_prediction <- predict_inactivation(t,
                         primary_model,
                         environment = "dynamic",
                         sec,
                         env_conditions)

    ## Prepare the output

    out <- list(
      method = method,
      algorithm = "regression",
      data = fit_data,
      guess = guess,
      known = known,
      primary_model = primary_model_name,
      fit_results = my_fit,
      best_prediction = best_prediction,
      sec_models = sec,
      env_conditions = env_conditions,
      niter = niter
      # logbase_logN = NULL,
      # approach_logN0 = NULL
    )

    class(out) <- c("InactivationFit", class(out))

    ## Return

    out

  } else if (method == "global") {

    ## TODO

  } else {
    stop("Unknown fitting method: ", method)
  }

}

#'
#' @export
#'
make_fitting_guess <- function(method) {

  if (method == "primary") {

  } else if (method == "two-step") {

  } else if (method == "one-step") {

  } else if (method == "dynamic") {

  } else if (method == "global") {

  } else {
    stop("Unknown fitting method: ", method)
  }

}


