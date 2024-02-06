
#' Residuals for fitting secondary models
#' 
#' The function is defined to be called within modFit or modMCMC
#' 
#' @param this_p named numeric vector of parameter values
#' @param model_name identifier of the secondary model as per [secondary_model_data()]
#' @param fit_data a tibble (or data.frame) defining the data as in [fit_inactivation_secondary()]
#' @param known a named numeric vector of known model parameters
#' @param output output format. Either 'vector' (default) or 'loglik'
#' 
#' @importFrom purrr imap_dfc
#' 
#' @returns a numeric vector of model residuals
#'
secondary_residuals <- function(this_p,
                                model_name,
                                fit_data,
                                known,
                                output = "vector"
                                ) {

  ## Put the parameters together

  p <- c(this_p, known)

  ## Make log-transformations in the parameters of the secondary model

  pars_log <- p[str_detect(names(p), "_log")]
  pars_log <- 10^pars_log
  names(pars_log) <- str_replace(names(pars_log), "log", "")
  p <- c(p, pars_log)
  p <- p[!str_detect(names(p), "_log")]  # remove the logs to avoid problems

  ## Some parameters are defined directly in "log" (without the "-"; e.g. logref)

  pars_log <- p[str_detect(names(p), "log")]
  pars_log <- 10^pars_log
  names(pars_log) <- str_replace(names(pars_log), "log", "")
  p <- c(p, pars_log)
  p <- p[!str_detect(names(p), "log")]  # remove the logs to avoid problems

  ## Calculate the effect of the environmental factors

  effects <- fit_data %>%
    select(-"my_par") %>%
    imap_dfc(
      ~ switch(model_name,  # TODO: I should consider putting this into a function (code repeated twice)
               Bigelow = sec_Bigelow(
                 x = .x,
                 xref = p[[paste0(.y, "_", "xref")]], z = p[[paste0(.y, "_", "z")]]
               ),
               genBigelow = sec_genBigelow(
                 x = .x,
                 xref = p[[paste0(.y, "_", "xref")]], z = p[[paste0(.y, "_", "z")]], n = p[[paste0(.y, "_", "n")]]
               ),
               Lineal = sec_Lineal(
                 x = .x,
                 xref = p[[paste0(.y, "_", "xref")]], b = p[[paste0(.y, "_", "b")]]
               ),
               Arhenius = sec_Arhenius(
                 x = .x,
                 xref = p[[paste0(.y, "_", "xref")]], Ea = p[[paste0(.y, "_", "Ea")]]
               ),
               LogExponential = sec_logExponential(
                 x = .x,
                 k = p[[paste0(.y, "_", "k")]], xc = p[[paste0(.y, "_", "Xcrit")]]
               ),
               stop(paste("Unknown model:", this_model))
      )
    )

  ## Put together with the reference value

  if (model_name %in% c("Bigelow", "genBigelow")) {

    pred <- 10^( log10(p[["ref"]]) + rowSums(effects) )

  } else {

    pred <- p[["ref"]] + rowSums(effects)

  }

  ## Output the residuals
  
  res <- fit_data$my_par - pred
  
  if (output == "loglik") {
    
    sum(res^2)
    
  } else {
    
    res
    
  }

  # # fit_data$my_par - pred
  # 
  # res <- fit_data$my_par - pred
  # sum(res^2)

}

#' Fitting secondary inactivation models
#' 
#' The function can fit secondary models with an arbitrary number of environmental factors
#' to a datset containing the values of some parameter observed in independent experiments
#' at constant environmental conditions. 
#' 
#' @details
#' 
#' The model variables are defined using the `formula` attribute. It must be a two-sided formula
#' with the left hand side defining a unique output variable (i.e., the parameter of the primary model) 
#' and the right hand side an arbitrary number of environmental factors (separated by `+`). 
#' For instance, `D ~ temp` would define a model for `D` as a function of `temp`, whereas
#' `D ~ temp + pH` would define a model for `D` as a function of both `temp` and `pH`. The name 
#' of these variables must match the column names in `fit_data`.
#' 
#' Each model parameter must be assigned either an initial guess (`guess`) or a fixed
#' value (`known`) as a named numeric vector. The conventions for the names are 
#' `factor-name` + `secondary-parameter`. For instance, the name `temperature_z` would provide
#' the guess for the z-value with respect to the factor
#' named "temperature". The names of the parameters must match those returned by
#' [secondary_model_data()] and the factor name must be included in `formula`. Note that the fitting can 
#' also be done in log-scale be including "log" before the parameter name: 
#' e.g., `temperature_logz`.
#' 
#' 
#' @param fit_data tibble (or data.frame) with the data for the fitting. Must have one column
#' with the parameter values of the primary model and as many columns as needed describing
#' environmental factors.
#' @param model_name A model identifier according to [secondary_model_data()]
#' @param formula A two-sided formula describing the output and input variables of the model. See details.
#' @param guess a named numeric vector with initial guesses for the model parameters. See details.
#' @param known a named numeric vector of parameter values that are considered known.
#' @param upper a named numeric vector with upper limits for the parameter estimates
#' @param lower a named numeric vector with lower limits for the parameter estimates
#' @param algorithm the fitting algorithm to use. Either 'regression' or 'MCMC'.
#' @param niter number of iterations of the MCMC algorithm. Ignored when algorithm = 'regression'
#' @param ... additional arguments for [modFit()] or [modMCMC()]
#'
#' @importFrom formula.tools lhs rhs get.vars
#' @importFrom FME modFit modMCMC
#'
#' @export
#'
fit_inactivation_secondary <- function(fit_data,
                                       model_name,
                                       formula = my_par ~ temp,
                                       guess,
                                       known,
                                       upper = NULL,
                                       lower = NULL,
                                       algorithm = "regression",
                                       niter = NULL,
                                       ...
                                       # check = NULL  # TODO 
                                       ) {

  ## Apply the formula

  y_col <- lhs(formula)

  vars <- all.vars(formula)
  vars <- vars[vars != y_col]

  fit_data <- select(fit_data,
                     my_par = y_col,
                     matches(vars)
  )
  
  ## Check the model parameters
  
  check_secondary_pars(model_name, c(guess, known), vars)
  
  ## Set up the bounds
  
  if (is.null(upper)) upper <- Inf
  if (is.null(lower)) lower <- -Inf
  
  ## Fit the model

  if (algorithm == "regression") {

    my_fit <- modFit(secondary_residuals,
                     unlist(guess),
                     fit_data = fit_data,
                     model_name = model_name,
                     known = unlist(known),
                     upper = upper,
                     lower = lower,
                     ...
    )

    ## Prepare the output

    out <- list(
      algorithm = "regression",
      data = fit_data,
      guess = guess,
      known = known,
      model_name = model_name,
      fit_results = my_fit,
      # best_prediction = best_prediction,  # TODO: implement this
      niter = NULL,
      par_name = y_col
    )

    class(out) <- c("SecondaryFit", class(out))

    ## Return

    out


  } else if (algorithm == "MCMC") {
    
    my_fit <- modMCMC(secondary_residuals,
                      unlist(guess),
                      fit_data = fit_data,
                      model_name = model_name,
                      known = unlist(known),
                      upper = upper,
                      lower = lower,
                      niter = niter,
                      output = "loglik",
                      ...
                      )
    
    ## Prepare the output
    
    out <- list(
      algorithm = "MCMC",
      data = fit_data,
      guess = guess,
      known = known,
      model_name = model_name,
      fit_results = my_fit,
      # best_prediction = best_prediction,  # TODO: implement this
      niter = niter,
      par_name = y_col
    )
    
    class(out) <- c("SecondaryFit", class(out))
    
    ## Return
    
    out

  } else {
    stop("Algorithm must be 'regression' or 'MCMC', got: ", algorithm)
  }

}

