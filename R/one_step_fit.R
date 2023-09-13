
#' Residuals for one-step fitting
#'
#' @importFrom stringr str_detect
#' @importFrom FME modCost
#' @importFrom dplyr pick
#'
#'
onestep_residuals <- function(this_p,
                              fit_data,
                              primary_model_name,
                              sec_models,
                              known
                              # approach_logN0
                              ) {

  ## Put the parameters together

  p <- c(this_p, known)

  ## Make log-transformations in the parameters of the secondary model

  pars_log <- p[str_detect(names(p), "_log")]
  pars_log <- 10^pars_log
  names(pars_log) <- str_replace(names(pars_log), "log", "")
  p <- c(p, pars_log)
  p <- p[!str_detect(names(p), "_log")]  # remove the logs to avoid problems

  ## Some parameters are defined directly in "log" (without the "-"; e.g. logN0)

  pars_log <- p[str_detect(names(p), "log")]
  pars_log <- 10^pars_log
  names(pars_log) <- str_replace(names(pars_log), "log", "")
  p <- c(p, pars_log)
  p <- p[!str_detect(names(p), "log")]  # remove the logs to avoid problems

  ## Prepare the secondary models

  sec <- convert_dynamic_guess(sec_models, p, c())

  cond <- select(fit_data, -time, -logN)

  primary_pars <- lapply(1:nrow(cond), function(i) {
    apply_secondary_models(cond[i,], sec) %>% as_tibble()
  }) %>%
    bind_rows() %>%
    mutate(time = fit_data$time)


  ## Calculate the logN  ----- I should put this into an independent function (same code twice)


  logN <- switch(primary_model_name,
                 Bigelow = iso_Bigelow(primary_pars$time,
                                       p[["N0"]],
                                       primary_pars$D),
                 Peleg = iso_Peleg(primary_pars$time,
                                   p[["N0"]],
                                   primary_pars$b,
                                   primary_pars$n),
                 Mafart = iso_Mafart(primary_pars$time,
                                     p[["N0"]],
                                     primary_pars$delta,
                                     primary_pars$p
                 ),
                 Geeraerd = iso_Geeraerd(primary_pars$time,
                                         p[["N0"]],
                                         primary_pars$D,
                                         primary_pars$Nres,
                                         primary_pars$SL
                 ),
                 Metselaar = iso_Metselaar(primary_pars$time,
                                           p[["Delta"]],
                                           p[["N0"]],
                                           primary_pars$D,
                                           primary_pars$p),
                 Weibull_2phase = iso_Weibull_2phase(primary_pars$time,
                                                     p[["N0"]],
                                                     primary_pars$delta1,
                                                     primary_pars$delta2,
                                                     primary_pars$p1,
                                                     primary_pars$p2,
                                                     primary_pars$alpha),
                 Trilinear = iso_Trilinear(primary_pars$time,
                                           p[["N0"]],
                                           primary_pars$SL,
                                           primary_pars$D,
                                           primary_pars$Nres),
                 Geeraerd_k = iso_Geeraerd_k(primary_pars$time,
                                             p[["N0"]],
                                             primary_pars$SL,
                                             primary_pars$k,
                                             primary_pars$Nres),
                 Geeraerd_k_noTail = iso_Geeraerd_noTail_k(primary_pars$time,
                                                           p[["N0"]],
                                                           primary_pars$SL,
                                                           primary_pars$k),
                 Geeraerd_noTail = iso_Geeraerd_noTail(primary_pars$time,
                                                       p[["N0"]],
                                                       primary_pars$SL,
                                                       primary_pars$D),
                 Geeraerd_k_noShoulder = iso_Geeraerd_noShoulder_k(primary_pars$time,
                                                                   p[["N0"]],
                                                                   primary_pars$k,
                                                                   primary_pars$Nres),
                 Geeraerd_noShoulder = iso_Geeraerd_noShoulder(primary_pars$time,
                                                               p[["N0"]],
                                                               primary_pars$D,
                                                               primary_pars$Nres),
                 stop(paste("Unknown model:", model_name))
  )

  # ## Calculate residuals
  # 
  # res <- fit_data$logN - logN
  # 
  # ## Return
  # 
  # res
  
  ## Calculate the cost
  
  modCost(
    model = data.frame(time = fit_data$time,
                       logN = logN),
    obs = select(fit_data, time, logN) %>% as.data.frame()
  )

}



#' One-step fitting of isothermal inactivation data
#' 
#' @param fit_data a tibble (or data.frame) with the data for the fitting.
#' @param model_name a model identifier according to [primary_model_data()]
#' @param start a named numeric vector of initial guesses for the parameters as
#' defined in [fit_inactivation()]
#' @param known a named numeric vector of known parameters
#' @param upper a named numeric vectors of upper bounds for the parameters
#' @param lower a named numeric vector of lower bounds for the parameters
#' @param secondary_models a nested list defining the secondary models for each
#' parameter as defined in [fit_inactivation()]
#' @param algorithm one of `"regression"` (default) or `"MCMC"`
#' @param niter number of MC iterations. Ignored if `algorithm == "regression`
#' @param ... additional arguments for [modFit()] or [modMCMC()]
#' 
#' @returns An instance of [InactivationFit]
#' 
#'
fit_onestep <- function(fit_data,
                        model_name,
                        start,
                        known,
                        upper = NULL,
                        lower = NULL,
                        secondary_models,
                        algorithm,
                        niter = NULL,
                        # approach_logN0
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

  # ## Deal with N0
  #
  # if (approach_logN0 == "unique") {
  #
  # } else if (approach_logN0 == "logS") {
  #
  #   my_cols <- fit_data %>%
  #     select(-time, -logN) %>%
  #     names()
  #
  #   fit_data <- fit_data %>%
  #     group_by(pick(my_cols)) %>%
  #     mutate(logN0 = mean(
  #       ifelse(time == 0, logN, NA),
  #       na.rm = TRUE
  #     )) %>%
  #     mutate(logS = logN - logN0) %>%
  #     ungroup()
  #
  # } else if (approach_logN0 == "different") {
  #
  # } else {
  #   stop(paste0("approach_logN0 must be 'unique', 'logS' or 'different', got: ", approach_logN0))
  # }
  
  ## Set up the bounds
  
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

  ## Fit the model

  if (algorithm == "regression") {

    my_fit <- modFit(onestep_residuals,
                     unlist(start),
                     fit_data = fit_data,
                     primary_model_name = model_name,
                     sec_models = secondary_models,
                     known = unlist(known),
                     upper = upper,
                     lower = lower,
                     # approach_logN0 = approach_logN0
                     ...
                     )
    
    ## Extract the secondary and primary models
    
    p <- c(coef(my_fit), unlist(known))
    
    primary_model <- list(model = model_name)
    
    initial <- p[str_detect(names(p), "N0")]
    aa <- initial
    names(aa) <- NULL
    primary_model[[names(initial)]] <- aa
    
    if (str_detect(model_name, "Geeraerd")) {
      initial <- p[str_detect(names(p), "C0")]
      aa <- initial
      names(aa) <- NULL
      primary_model[[names(initial)]] <- aa
    }
    
    sec <- convert_dynamic_guess(secondary_models, p, c())
    
    ## Prepare the output
    
    out <- list(
      approach = "one-step",
      algorithm = "regression",
      data = fit_data,
      guess = start,
      known = known,
      primary_model = model_name,
      fit_results = my_fit,
      best_prediction = NA,
      sec_models = sec,
      env_conditions = NULL,
      niter = NULL
      # approach_logN0 = NULL
    )
    
    class(out) <- c("InactivationFit", class(out))


  } else if (algorithm == "MCMC") {
    
    my_fit <- modMCMC(onestep_residuals,
                     unlist(start),
                     fit_data = fit_data,
                     primary_model_name = model_name,
                     sec_models = secondary_models,
                     known = unlist(known),
                     upper = upper,
                     lower = lower,
                     niter = niter,
                     # approach_logN0 = approach_logN0
                     ...
    )
    
    ## Extract the secondary and primary models
    
    p <- c(my_fit$bestpar, unlist(known))
    
    primary_model <- list(model = model_name)
    
    initial <- p[str_detect(names(p), "N0")]
    aa <- initial
    names(aa) <- NULL
    primary_model[[names(initial)]] <- aa
    
    if (str_detect(model_name, "Geeraerd")) {
      initial <- p[str_detect(names(p), "C0")]
      aa <- initial
      names(aa) <- NULL
      primary_model[[names(initial)]] <- aa
    }
    
    sec <- convert_dynamic_guess(secondary_models, p, c())
    
    ## Prepare the output
    
    out <- list(
      approach = "one-step",
      algorithm = "MCMC",
      data = fit_data,
      guess = start,
      known = known,
      primary_model = model_name,
      fit_results = my_fit,
      best_prediction = NA,
      sec_models = sec,
      env_conditions = NULL,
      niter = niter
      # approach_logN0 = NULL
    )
    
    class(out) <- c("InactivationFit", class(out))
    

  } else {
    stop("Algorithm must be 'regression' or 'MCMC', got: ", algorithm)
  }

  ## Output

  out

}
