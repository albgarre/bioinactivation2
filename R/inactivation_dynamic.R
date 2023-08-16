

#'
predict_dynamic_inactivation <- function(times,
                                         primary_model,
                                         env_conditions,
                                         secondary_models,
                                         ...,
                                         # check = TRUE,
                                         logbase_logN = 10
                                         # logbase_mu = logbase_logN,
                                         # formula = . ~ time
) {

  # ## Apply the formula
  #
  # x_col <- rhs(formula)
  #
  # env_conditions <- rename(env_conditions, time = x_col)
  #
  # ## Check model parameters
  #
  # if (isTRUE(check)) {
  #
  #   sec_model_names <- secondary_models %>%
  #     map(~ .$model) %>%
  #     unlist()
  #
  #   check_pars_names <- lapply(names(secondary_models), function(each_factor) {
  #
  #     this_model <- secondary_models[[each_factor]]
  #     this_model$model <- NULL
  #     paste0(each_factor, "_", names(this_model))
  #   }) %>%
  #     unlist()
  #
  #   check_pars <- rep(1, length(check_pars_names))  # I do not check values, so it does not matter
  #   names(check_pars) <- check_pars_names
  #
  #
  #   check_secondary_pars(check_pars, unlist(primary_pars), sec_model_names,
  #                        primary_pars = c("mu_opt", "N0", "Nmax", "Q0"))
  #
  # }

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
             env_interpolator = my_env, secondary_models = secondary_models) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(logN = log(.data$N, base = logbase_logN))

  ## Return

  out

}




