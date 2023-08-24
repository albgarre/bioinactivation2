
#' Showing initial guesses for primary models
#' 
#' @param fit_data data for the fit. A tibble with two columns: 'time' and 'logN'
#' @param primary_model_name model identifier as per [primary_model_data()]
#' @param guess a named numeric vector with the initial guess, as in [predict_inactivation()]
#' 
#' @returns an instance of ggplot comparing the prediction against the data.
#' 
show_guess_primary <- function(fit_data, primary_model_name, guess) {
  
  
  ## Calculate the prediction
  
  primary_model <- as.list(guess)
  primary_model$model <- primary_model_name
  t <- seq(0, max(fit_data$time), length = 1000)
  
  pred <- predict_inactivation(t, primary_model)
  
  ## Make the plot
  
  plot(pred) +
    geom_point(aes(x = time, y = logN), data = fit_data)
  
}

#' Showing initial guesses for dynamic conditions
#' 
#' @inheritParams show_guess_primary
#' @param sec_models a nested list defining the secondary models as in [fit_inactivation()]
#' @param env_conditions a tibble (or data.frame) describing the environmental conditions as
#' in [fit_inactivation()]
#' 
#' @importFrom stringr str_detect
#' 
#' @returns an instance of ggplot comparing the prediction against the data
#' 
show_guess_dynamic <- function(fit_data, primary_model_name, guess,
                               sec_models, env_conditions) {
  
  ## Calculate the prediction
  
  primary_model <- list(model = primary_model_name)
  
  initial <- guess[str_detect(names(guess), "N0")]
  aa <- initial
  names(aa) <- NULL
  primary_model[[names(initial)]] <- aa
  
  if (str_detect(primary_model_name, "Geeraerd")) {
    initial <- guess[str_detect(names(guess), "C0")]
    aa <- initial
    names(aa) <- NULL
    primary_model[[names(initial)]] <- aa
  }
  
  t <- seq(0, max(fit_data$time), length = 1000)
  
  sec <- convert_dynamic_guess(sec_models, guess, c())
  
  pred <- predict_inactivation(t,
                               primary_model,
                               environment = "dynamic",
                               sec,
                               env_conditions)
  
  ## Make the plot
  
  plot(pred) +
    geom_point(aes(x = time, y = logN), data = fit_data)
  
}

#' Checking guesses for inactivation models
#' 
#' Makes a plot comparing the prediction corresponding to some parameter values
#' against the experimental data
#' 
#' @param method the fitting method. One of 'primary', 'two-steps', 'one-step', dynamic'
#' or 'global'
#' @param fit_data the data used for the fit, as defined in [fit_inactivation()]
#' @param primary_model_name a model identifier as in [fit_inactivation()]
#' @param guess a named numeric vector defining the model parameters as in [fit_inactivation()]
#' @param secondary_models a nested list defining the secondary models for each parameter
#' as in [fit_inactivation()]
#' @param env_conditions a tibble (or data.frame) describing the environmental conditions
#' during the treatment as in [fit_inactivation()]
#' @param formula a formula defining the variables of the primary model. By default, logN ~ time 
#' 
#' @returns an instance of ggplot comparing the prediction against the data.
#' 
#' @export
#' 
check_inactivation_guess <- function(method,
                                     fit_data,
                                     primary_model_name,
                                     guess,
                                     # approach_logN0 = "unique",  # or "logS" or "different",
                                     secondary_models = NULL,
                                     env_conditions = NULL,
                                     formula = logN ~ time
                                     ) {
  
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
  
  if (method == "primary") {
    
    show_guess_primary(fit_data, primary_model_name, guess)
    
  } else if (method == "two-steps") {
    
    ## AA
    
    
  } else if (method == "one-step") {
    
    ## AA
    
    
  } else if (method == "dynamic") {
    
    show_guess_dynamic(fit_data, primary_model_name, guess,
                       sec_models, env_conditions)
    
  } else if (method == "global") {
    
    ## TODO
    
  } else {
    stop("Unknown fitting method: ", method)
  }
  
}

#' Checking guesses for secondary inactivation models
#' 
#' @param fit_data a tibble (or data.frame) with the data to use for the fit, as in
#' [fit_inactivation_secondary()]
#' @param model_name a secondary model name as in [secondary_model_data()]
#' @param guess a named numeric vector with the model parameters
#' @param formula a formula defining the model variables as in [fit_inactivation_secondary()]
#' 
#' @importFrom cowplot plot_grid
#' 
#' @returns an instance of ggplot with two subplots: an observed vs predicted plot 
#' and a predicted vs residuals plot
#' 
#' @export
#' 
check_secondary_guess <- function(fit_data,
                                  model_name,
                                  guess,
                                  formula = my_par ~ temp
                                  ) {
  
  ## Apply the formula
  
  y_col <- lhs(formula)
  
  vars <- all.vars(formula)
  vars <- vars[vars != y_col]
  
  fit_data <- select(fit_data,
                     my_par = y_col,
                     matches(vars)
  )
  
  # ## Check the model parameters
  # 
  # check_secondary_pars(model_name, c(start, known), vars)
  
  ## Calculate the residuals
  
  res <- secondary_residuals(guess, model_name, fit_data, c())
  
  fit_data <- fit_data %>%
    mutate(res = res,
           pred = my_par - res)
  
  p1 <- ggplot(fit_data, aes(x = my_par, y = pred)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = 2)
  
  p2 <- ggplot(fit_data, aes(x = pred, y = res)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    geom_hline(yintercept = 0, linetype = 2) 
  
  plot_grid(p1, p2)
  
}























