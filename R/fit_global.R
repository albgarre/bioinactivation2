

#' Residuals of multiple dynamic predictions
#' 
#' @param experiment_data a nested list with the experimental data. Each entry describes
#' one experiment as a list with two elements: data and conditions. `data` is a tibble
#' with two columns: time and logN. `conditions` is a tibble with one column named time
#' and as many additional columns as environmental factors.
#'
#' @return an instance of `modCost.
#'
#' @importFrom FME modCost
#'
get_multi_dyna_residuals <- function(this_p, 
                                     experiment_data,
                                     known, 
                                     primary_model_name,
                                     sec_models) {
  
  old_cost <- NULL
  
  for (each_experiment in experiment_data) {
    
    my_cost <- dynamic_residuals(this_p = this_p,
                                 fit_data = each_experiment$data,
                                 primary_model_name = primary_model_name,
                                 sec_models = sec_models,
                                 known = known,
                                 env_conditions = each_experiment$conditions)
    
    
    old_cost <- my_cost
    
  }
  
  my_cost
  
}


#' Fitting a single model to several (dynamic) experiments
#' 
fit_multiple_inactivation <- function(fit_data,
                                      primary_model_name,
                                      guess,
                                      known,
                                      upper = NULL,
                                      lower = NULL,
                                      secondary_models,
                                      algorithm = "regression",
                                      env_conditions,
                                      niter = NULL,
                                      ...
                                      # check = TRUE,
                                      # formula = logN ~ time
                                      ) {
  
  
  # ## Check the model parameters
  #
  # if (isTRUE(check)) {
  #
  #   check_primary_pars(model_name, c(starting_point, known_pars))
  #
  # }
  #
  # ## Apply the formula
  #
  # if (length(get.vars(formula)) > 2) {
  #   stop("Only formulas with 2 terms are supported.")
  # }
  #
  # y_col <- lhs(formula)
  # x_col <- rhs(formula)
  #
  # fit_data <- select(fit_data,
  #                    time = x_col,
  #                    logN = y_col
  # )
  
  ## Put the environmental and the microbial datas together
  
  if (is.null(names(fit_data))) {  # In case the list were unnamed
    names(fit_data) <- paste0("exp_", 1:length(fit_data))
    names(env_conditions) <- paste0("exp_", 1:length(fit_data))
  }
  
  my_data <- names(fit_data) %>%
    map(
      ~ list(data = fit_data[[.]],
             conditions = env_conditions[[.]]
             )
      )
  
  names(my_data) <- names(fit_data)
  
  ## Set up the bounds
  
  if (is.null(upper)) upper <- Inf
  if (is.null(lower)) lower <- -Inf
  
  ## Fit the model
  
  if (algorithm == "regression") {
    
    my_fit <- modFit(
      get_multi_dyna_residuals,
      guess, 
      experiment_data = my_data,
      known = known, 
      primary_model_name = primary_model_name,
      sec_models = sec_models,
      upper = upper,
      lower = lower,
      ...
      )
    
    ## Prepare the output
    
    out <- list(
      approach = "global",
      algorithm = "regression",
      data = my_data,
      guess = guess,
      known = known,
      primary_model = primary_model_name,
      fit_results = my_fit,
      # best_prediction = best_prediction,
      sec_models = sec_models,
      # env_conditions = env_conditions,
      niter = NULL
    )
    
    class(out) <- c("InactivationFit", class(out))
    
    
  } else if (algorithm == "MCMC") {
    
    my_fit <- modMCMC(
      get_multi_dyna_residuals,
      guess, 
      experiment_data = my_data,
      known = known, 
      primary_model_name = primary_model_name,
      sec_models = sec_models,
      upper = upper,
      lower = lower,
      niter = niter,
      ...
    )
    
    ## Prepare the output
    
    out <- list(
      approach = "global",
      algorithm = "MCMC",
      data = my_data,
      guess = guess,
      known = known,
      primary_model = primary_model_name,
      fit_results = my_fit,
      # best_prediction = best_prediction,
      sec_models = sec_models,
      # env_conditions = env_conditions,
      niter = niter
    )
    
    class(out) <- c("InactivationFit", class(out))
    
  } else {
    stop("Algorithm must be 'regression' or 'MCMC', got: ", algorithm)
  }
  
  ## Output
  
  out
  
}






