
#'
#'
secondary_residuals <- function(this_p,
                                model_name,
                                fit_data,
                                # sec_models,
                                known
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

  # ## Prepare the secondary models
  #
  # sec_models
  #
  # sec <- convert_dynamic_guess(sec_models, p, c())
  #
  # cond <- select(fit_data, -x)
  #
  # primary_pars <- lapply(1:nrow(cond), function(i) {
  #   apply_secondary_models(cond[i,], sec) %>% as_tibble()
  # }) %>%
  #   bind_rows() %>%
  #   mutate(time = fit_data$time)

  ## Calculate the effect of the environmental factors

  effects <- fit_data %>%
    select(-my_par) %>%
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

  fit_data$my_par - pred

}

#' Fits a secondary inactivation model using different methods from predictive microbiology
#'
#' @importFrom formula.tools lhs rhs get.vars
#'
#' @export
#'
fit_inactivation_secondary <- function(fit_data,
                                       model_name,
                                       start,
                                       known,
                                       upper = NULL,
                                       lower = NULL,
                                       algorithm = "regression",
                                       niter = NULL,
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
  
  ## Check the model parameters
  
  check_secondary_pars(model_name, c(start, known), vars)
  
  ## Fit the model

  if (algorithm == "regression") {

    my_fit <- modFit(secondary_residuals,
                     unlist(start),
                     fit_data = fit_data,
                     model_name = model_name,
                     # sec_models = secondary_models,
                     known = unlist(known)
                     # approach_logN0 = approach_logN0
                     # logbase_logN = logbase_logN,
                     # ...
    )

    my_fit

    ## Prepare the output

    out <- list(
      algorithm = "regression",
      data = fit_data,
      guess = guess,
      known = known,
      model_name = model_name,
      fit_results = my_fit,
      # best_prediction = best_prediction,  # TODO: implement this
      # sec_models = secondary_models,
      niter = NULL,
      logbase_logN = NULL,
      par_name = y_col
    )

    class(out) <- c("SecondaryFit", class(out))

    ## Return

    out


  } else if (algorithm == "MCMC") {

  } else {
    stop("Algorithm must be 'regression' or 'MCMC', got: ", algorithm)
  }




}

