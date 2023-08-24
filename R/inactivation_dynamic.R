

#' Predictions under dynamic environmental conditions
#' 
#' The solution is calculated by solving the differential equation using numerical
#' methods through the ode function from package deSolve.
#' 
#' @param times a numeric vector of time points where to calculate the solution
#' @param primary_model a list defining the primary model as in [predict_inactivation()]
#' @param env_conditions a tibble (or data.frame) describing the variation of the environmental
#' conditions as in [predict_inactivation()]
#' @param secondary_models a nested list describing the secondary models for each parameter,
#' as defined in [predict_inactivation]
#' @param ... additional arguments passed to ode
#' @param check whether to do some basic model checks. `TRUE` by default
#' 
#' @importFrom deSolve ode
#' 
#' 
#' 
predict_dynamic_inactivation <- function(times,
                                         primary_model,
                                         env_conditions,
                                         secondary_models,
                                         ...,
                                         check = TRUE
                                         # formula = . ~ time
                                         ) {

  # ## Apply the formula
  #
  # x_col <- rhs(formula)
  #
  # env_conditions <- rename(env_conditions, time = x_col)
  #
  ## Check model parameters

  if (isTRUE(check)) {
    
    check_dynamic_pars(primary_model, secondary_models)

  }
  
  ## Approximate the env conditions

  my_env <- approx_env(env_conditions)

  ## Prepare the vector of initial conditions

  yini <- c(N = primary_model$N0)

  if ("C0" %in% names(primary_model)) {
    yini[["C"]] <- primary_model$C0
  }

  ## Assign names to the secondary models

  par_names <- secondary_models %>% map_chr(~ .$par)

  secondary_models <- secondary_models %>%
    set_names(par_names)

  ## Pick the model

  ode_model <- switch(primary_model$model,
                      Bigelow = dyna_Bigelow,
                      Mafart = dyna_Mafart,
                      Peleg = dyna_Peleg,
                      Geeraerd = dyna_Geeraerd,
                      Geeraerd_noTail = dyna_Geeraerd_noTail,
                      Geeraerd_noShoulder = dyna_Geeraerd_noSL,
                      Geeraerd_k = dyna_Geeraerd_k,
                      Geeraerd_k_noTail = dyna_Geeraerd_noTail_k,
                      Geeraerd_k_noShoulder = dyna_Geeraerd_noSL_k,
                      
                      stop("Uknown model: ", primary_model$model)
                      )

  ## Call ode

  out <- ode(yini, times, ode_model, primary_model,
             env_interpolator = my_env, 
             secondary_models = secondary_models,
             ...
             ) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(logN = log10(.data$N))
    # mutate(logN = log(.data$N, base = logbase_logN))

  ## Return

  out

}




