
#' Residuals for dynamic fitting
#'
#' @importFrom stringr str_detect
#' @importFrom FME modCost
#'
#'
dynamic_residuals <- function(this_p,
                              fit_data,
                              primary_model_name,
                              sec_models,
                              known
                              ) {

  # browser()

  ## Put the parameters together

  p <- c(this_p, known)

  ## Prepare the primary model

  primary_model <- list(model = primary_model_name)

  ## Extract the initial count

  initial <- p[str_detect(names(p), "N0")]
  aa <- initial
  names(aa) <- NULL
  primary_model[[names(initial)]] <- aa

  ## If Geeraerd, extract C0

  if (str_detect(primary_model_name, "Geeraerd")) {
    initial <- p[str_detect(names(p), "C0")]
    aa <- initial
    names(aa) <- NULL
    primary_model[[names(initial)]] <- aa
  }

  ## Prepare the secondary models

  sec <- convert_dynamic_guess(sec_models, this_p, known)

  ## Make the prediction

  pred <- predict_inactivation(fit_data$time,
                               primary_model,
                               environment = "dynamic",
                               sec,
                               env_conditions)

  ## Calculate residuals

  pred$simulation %>%
    select(time, logN) %>%
    as.data.frame() %>%
    modCost(
      model = .,
      obs = as.data.frame(fit_data)
    )
}

#' AA
#'
#' @importFrom FME modFit modCost
#'
fit_dynamic_inactivation <- function(fit_data,
                                     model_name,
                                     start,
                                     known,
                                     upper,
                                     lower,
                                     secondary_models,
                                     algorithm,
                                     env_conditions,
                                     niter
                                     # ...,
                                     # check = TRUE,
                                     # logbase_logN = 10,  # TODO
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

  ## Fit the model

  # browser()

  if (algorithm == "regression") {

    my_fit <- modFit(dynamic_residuals,
                     unlist(start),
                     fit_data = fit_data,
                     primary_model_name = model_name,
                     sec_models = secondary_models,
                     known = unlist(known)
                     # logbase_logN = logbase_logN,
                     # ...
    )


  } else if (algorithm == "MCMC") {

  } else {
    stop("Algorithm must be 'regression' or 'MCMC', got: ", algorithm)
  }

  ## Output

  my_fit

}
