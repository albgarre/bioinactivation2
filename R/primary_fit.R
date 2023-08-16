
#'
#'
primary_residuals <- function(this_p,
                              fit_data,
                              model_name,
                              known
                              # logbase_logN = 10
) {

  ## Make the prediction

  times <- sort(unique(fit_data$time))

  pars <- c(this_p, known)
  my_model <- as.list(pars)
  my_model$model <- model_name

  pred <- predict_inactivation(times,
                       my_model,
                       check = FALSE
                       # logbase_mu = logbase_mu,
                       # logbase_logN = logbase_logN
                       )

  ## TODO: change the pred when I define the class

  ## Return the residuals

  modCost(model = as.data.frame(pred$simulation),
          obs = as.data.frame(fit_data)
          )

}


#' AA
#'
#' @importFrom FME modFit modCost
#'
fit_primary <- function(fit_data,
                        model_name,
                        start,
                        known,
                        upper = NULL,
                        lower = NULL,
                        # approach_logN0 = "unique",  # or "logS" or "different",
                        # secondary_models = NULL
                        # algorithm = "regression",
                        # env_conditions = NULL,
                        # niter = NULL,
                        # ...,
                        check = TRUE
                        # logbase_mu = logbase_logN,
                        # logbase_logN = 10,  # TODO
                        # formula = logN ~ time
) {

  ## Check the model parameters

  if (isTRUE(check)) {

    check_primary_pars(model_name, c(start, known))

  }
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

  my_fit <- modFit(primary_residuals,
                   unlist(start),
                   fit_data = fit_data,
                   model_name = model_name,
                   known = unlist(known)
                   # logbase_logN = logbase_logN,
                   # ...
                   )

  ## Output

  my_fit

}


