
#' Fitting inactivation models
#' 
#' This is a top-level function to fit inactivation models using different approachs 
#' from predictive microbiology. The function can fit primary models to a set of data
#' gathered (in principle) under constant environmental conditions, a unique model
#' (primary + secondary model) to a set of independent experiments at different 
#' (constant) environmental conditions, a unique model (primary + secondary model)
#' to data obtained under dynamic environmental conditions or a unique model (
#' primary + secondary model) to a set of independent experiments under dynamic (or not)
#' environmental conditions. Please see the sections below for additional information
#' on how each fitting must be defined.
#' 
#' @section Fitting primary models:
#' Setting `approach = "primary"` fits a primary inactivation model to a dataset describing
#' the inactivation of a population. In principles, this approach is valid to data
#' gathered under constant environmental conditions.
#' 
#' In this case, the data is defined through `fit_data` as a tibble
#' (or data.frame) with two columns: one describing the elapsed time (`time`) and a second
#' one describing the variation in the log-microbial concentration (`logN`). Note that
#' these column names can be changed using the `formula` argument.
#' 
#' The function allows fixing any model parameter before doing the fit using the `known`
#' argument. This is useful for cases where the parameter is not "biological" (e.g., `xref`)
#' or when there are  identifiability issues. Initial guesses must be provided to
#' the rest of the model parameters, due to the requirements of the fitting algorithms
#' implemented in the function. In either case, the parameters are defined as a
#' named numeric vector with names matching those returned by [primary_model_data()]. 
#' The fitting can also be done on the log-transformed parameter. To do this, you
#' just have to include "log" before the parameter name (e.g., `c(logD = 2)` instead
#' of `c(D = 2)`). Please see the examples below for more info.
#' Note that [make_guess_primary()] and [check_inactivation_guess()] can aid in the
#' definition of initial guesses for the parameters.
#' 
#' @section Defining secondary models for fitting:
#' Every fitting approach but `approach = "primary"` requires the definition of secondary
#' models. This is done in two steps: defining the model structure using `secondary_models`
#' and defining initial guesses for the parameter values using `guess` (or known values
#' using `known`).
#' 
#' The structure of the secondary models is defined using the `secondary_models` argument.
#' This is a nested list, where each element is a list that defines the secondary
#' model for one parameter of the primary model. These lists must have three elements:
#' * `par`: the name of the model parameter (as in [primary_model_data()])
#' * `model`: a model identifier (as per [secondary_model_data()])
#' * `depends_on`: a character vector with the names of the environmental factors that 
#' affect this parameter. Note that these names must match those in `fit_data`. If the
#' parameter does not depend on any environmental condition, an empty vector (`c()`) 
#' should be passed here. Please see the examples below for further information on
#' how to define secondary models.
#' 
#' Each model parameter must be assigned either an initial guess (`guess`) or a fixed
#' value (`known`) as a named numeric vector. The conventions for the names are `primary-parameter`
#' + `_` + `factor-name` + `secondary-parameter`. For instance, the name `delta_temperature_z` would provide
#' the guess for the z-value of the secondary model of parameter `delta` with respect to the factor
#' named "temperature". The names of the parameters must match those returned by
#' [primary_model_data()] and [secondary_model_data()]. Note that the fitting can 
#' also be done in log-scale be including "log" before the parameter name: 
#' e.g., `delta_temperature_logz`.
#' 
#' 
#' @section Fitting models to several experiments using the one-step approach: 
#' Setting `approach = "one-step"` fits a unique model (defined as a primary + secondary model)
#' to a dataset containing the inactivation of a population obtained under independent
#' experiments at (different) constant environmental conditions. In this approach,
#' the complete model (i.e., the combination of primary and secondary model) is
#' fitted to the complete dataset in a single step using nonlinear regression (or an
#' equivalent MCMC algorithm). 
#' 
#' The data is defined through `fit_data` as a tibble (or data.frame) with at least three columns.
#' The first two ones describe the elapsed time of the experiment (`time`) and the 
#' log-microbial concentration observed (`logN`). Note that these names can be changed 
#' using the `formula` argument. Then, additional columns define the values of the environmental
#' conditions of the experiment. In this sense, the function admits an arbitrary number
#' of environmental conditions. Note that the function does not make any consideration
#' regarding the number of experiments (it can be seen as each row being predicted independently),
#' so there is no requirement regarding the number of data points per experiment (
#' unlike for the two-steps fitting).
#' 
#' 
#' 
#' @section Fitting models to several experiments using the two-steps approach: 
#' Setting `approach = "two-steps"` fits a unique model (defined as a primary + secondary model)
#' to a dataset containing the inactivation of a population obtained under independent
#' experiments at (different) constant environmental conditions. In this approach,
#' a primary model is fitted to each experiment. Then, on a second step, a secondary
#' model is fitted to the parameters of the primary model estimated on the first step.
#' 
#' The data is defined through `fit_data` as a tibble (or data.frame) with the same conventions
#' as for the one-step fitting. That is, two columns  describe the elapsed time of the experiment (`time`) and the 
#' log-microbial concentration observed (`logN`). Note that these names can be changed 
#' using the `formula` argument. Then, additional columns define the values of the environmental
#' conditions of the experiment. In this sense, the function admits an arbitrary number
#' of environmental conditions, with the function identifying each experiment as rows
#' sharing the same environmental conditions. Hence, it is important that each experiment
#' includes enough data points to fit a primary model. Otherwise, the fitting algorithm is
#' unlikely to converge.
#' 
#' @section Fitting models to a dynamic experiment:
#' Setting `approach = "dynamic"` fits a unique model (a primary + secondary models)
#' to a dataset containing the inactivation of a population under an experiment with
#' varying environmental conditions. This is done by estimating the parameters of the model
#' defined directly as a differential equation (which is solved numerically).
#' 
#' The microbial data is defined through `fit_data` with the same conventions as for
#' fitting primary models: a tibble (or data.frame) with two columns describing the
#' treatment time (`time`) and the observed log-microbial concentration (`logN`). 
#' Then, the environmental conditions during the experiment are defined using the
#' `env_conditions` argument. This is a tibble (or data.frame) with one column named
#' `time` and as many additional columns as environmental conditions included in the model.
#' Model predictions are made by linear interpolation between the values provided in
#' `env_conditions`.
#' 
#' @section Fitting models to a set of dynamic experiments:
#' aa 
#'  
#' 
#' 
#' @param approach approach for model fitting. One of 'primary', 'two-steps', 'one-step',
#' 'dynamic' or 'global'.
#' @param fit_data a tibble with the data to use for model fitting. See sections below for further information.
#' @param primary_model_name a model identifier compatible with [primary_model_data()] or
#' [dynamic_model_data()] depending on the fitting approach.
#' @param guess a named numeric vector of initial guesses for the model parameters. See details
#' for conventions on how to define them
#' @param known a named numeric vector of parameter values that are considered known.
#' @param upper a named numeric vector with upper limits for the parameter estimates
#' @param lower a named numeric vector with lower limits for the parameter estimates
#' @param secondary_models a nested list defining the structure of the secondary models. 
#' See Details for further info.
#' @param algorithm the fitting algorithm to use. Either 'regression' or 'MCMC'.
#' @param env_conditions a tibble describing the change in the environmental conditions. Must have a column
#' named 'time' (although this can be changed using the formula argument) and as many
#' additional columns as environmental factors included in the model (with names matching
#' those in the `secondary_models` argument). 
#' @param niter number of iterations of the MCMC algorithm. Ignored when algorithm = 'regression'
#' @param check whether to do some basic checks of input consistency. By default, TRUE.
#' @param ... additional arguments passed to [modMCMC()] or [modFit()]
#' @param formula a formula describing the variables names for the primary model. By default, `logN ~ time` 
#' 
#' @returns an instance of [InactivationFit] with the fitted model
#' 
#' 
#' @export
#'
fit_inactivation <- function(approach,
                             fit_data,
                             primary_model_name,
                             guess = NULL,
                             known = c(),
                             upper = NULL,
                             lower = NULL,
                             # approach_logN0 = "unique",  # or "logS" or "different",
                             secondary_models = NULL,
                             algorithm = "regression",
                             env_conditions = NULL,
                             niter = NULL,
                             ...,
                             check = TRUE,
                             formula = logN ~ time
                             ) {

  if (approach == "primary") {
    
    ## Return an error if the initial guess is missing
    
    if ( is.null(guess) ) stop("Please provide an initial guess for the model parameters")
    
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

    ## Fit the model

    out <-  fit_primary(fit_data, 
                           primary_model_name, 
                           guess, known,
                           check = check,
                           upper = upper,
                           lower = lower,
                           algorithm = algorithm,
                           niter = niter,
                           ...
                           )

    ## Return

    out

  } else if (approach == "two-steps") {

    out <- fit_two_step(
      fit_data = fit_data,
      model_name = primary_model_name,
      known = known,
      secondary_models = secondary_models
      # approach_logN0 = approach_logN0
    )

    ## Return

    out



  } else if (approach == "one-step") {
    
    ## Return an error if the initial guess is missing
    
    if ( is.null(guess) ) stop("Please provide an initial guess for the model parameters")

    ## Fit the model

    out <- fit_onestep(
      fit_data = fit_data,
      model_name = primary_model_name,
      start = guess,
      known = known,
      upper = upper,
      lower = lower,
      secondary_models = secondary_models,
      algorithm = algorithm,
      niter = niter
      # approach_logN0 = approach_logN0
    )

    ## Return

    out


  } else if (approach == "dynamic") {
    
    ## Return an error if the initial guess is missing
    
    if ( is.null(guess) ) stop("Please provide an initial guess for the model parameters")

    ## Fit the model

    out <- fit_dynamic_inactivation(
      fit_data = fit_data,
      model_name = primary_model_name,
      start = guess,
      known = known,
      upper = upper,
      lower = lower,
      secondary_models = secondary_models,
      algorithm = algorithm,
      env_conditions = env_conditions,
      niter = niter
    )

    ## Return

    out

  } else if (approach == "global") {
    
    ## Return an error if the initial guess is missing
    
    if ( is.null(guess) ) stop("Please provide an initial guess for the model parameters")

    out <- fit_multiple_inactivation(fit_data,
                                        primary_model_name,
                                        guess,
                                        known,
                                        upper,
                                        lower,
                                        secondary_models,
                                        algorithm,
                                        env_conditions,
                                        niter)
    
    ## Return
    
    out

  } else {
    stop("Unknown fitting approach: ", approach)
  }

}

