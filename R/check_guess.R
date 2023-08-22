
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





