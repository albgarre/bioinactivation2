
#' Prediction of microbial inactivation
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' This function provides a top-level interface for predicting population inactivation. 
#' Predictions can be made either under constant or dynamic environmental conditions. 
#' See below for details on the calculations.
#' 
#' @details 
#' To ease data input, parameters can be defined with or without a log transformation.
#' For that, one has to include 'log' before the parameter name in parameter definition.
#' For instance, c(D=6) would define a D-value of 6 (units), whereas c(logD=6) would define
#' a D-value of 6 log-units. Please see examples and the package vignettes for further
#' examples.
#' 
#' @section Predictions in constant environments:
#' Predictions under constant environments are calculated using only the algebraic
#' form of the primary models (see vignette for details).
#' Consequently, the arguments "secondary_models" and "env_conditions" are ignored.
#' If these were passed, the function would return a warning.
#' 
#' The inactivation model is defined through the "primary_model" argument using a named list.
#' One of the list elements must be named "model" and must take take one of the valid 
#' keys returned by [primary_model_data()]. The remaining entries of the list define the
#' values of the parameters of the selected model. A list of valid keys can be retrieved
#' using [primary_model_data()] (see example below). Note that parameters can be defined
#' in different scales (see Details).
#' 
#' @section Predictions in dynamic environments:
#' Predictions under dynamic environments are calculated by solving numerically
#' the differential equation of the primary inactivation model. The effect of
#' changes in the environmental conditions in the inactivation rate are calculated
#' according to the secondary model defined by the user. Hence, predictions under
#' dynamic conditions require the definition of both primary and secondary models.
#' 
#' The dynamic environmental conditions are defined using a tibble (or data.frame)
#' through the "env_conditions" argument. It must include one column named "time" stating the elapsed
#' time and as many additional columns as environmental conditions included in the model
#' For values of time not included in the tibble, the values of the environmental conditions
#' are calculated by linear interpolation.
#' 
#' Primary models are defined as a named list through the "primary_model" argument. It is
#' defined in a similar way as for predictions under constant conditions. It is a named list
#' with an element named 'model' defining the primary model (according to `dynamic_model_data()`).
#' Additional arguments define the initial concentration (N0) and, for the Geeraerd model, the initial
#' value of the ideal substance C.
#' 
#' Secondary models are defined as a nested list through the "secondary_models" argument.
#' The list must have one entry per model parameter (according to `dynamic_model_data()`).
#' Then, the secondary model for each parameter is defined as a list with several arguments:
#' * par: the parameter defined (according to `dynamic_model_data()`).
#' * model: the secondary model used for this parameter (according to `secondary_model_data()`)
#' * ref: the value of the parameter for reference conditions. Note that defining the parameter
#' name as 'log'+name.
#' * As many additional entries as environmental factors that affect the parameter. Each
#' of these entries must be a named numeric vector defining the model parameters as
#' per `secondary_model_data()`. For instance, pH = c("xref" = 7, z = 2) would define
#' a dependency with respect to pH with x_ref of 7 and z-value of 2. Please see the
#' examples and the vignettes for additional data.
#' 
#' @param times numeric vector of time points for making the predictions
#' @param primary_model  named list defining the values of the parameters of the primary inactivation model
#' @param environment type of environment. Either "constant" (default) or "dynamic" (see below for details 
#' on the calculations for each condition)
#' @param secondary_models a nested list describing the secondary models. See below for details
#' @param env_conditions Tibble describing the variation of the environmental
#' conditions for dynamic experiments. It must have with the elapsed time (named `time` 
#' by default; can be changed with the "formula" argument), 
#' and as many additional columns as environmental factors. Ignored for "constant" environments.
#' @param ... Additional arguments for [ode()].
#' @param check Whether to check the validity of the models. `TRUE` by default.
#' @param formula An object of class "formula" describing the x variable for predictions 
#' under dynamic conditions. `. ~ time` as a default.
#' 
#' @return An instance of [InactivationPrediction].
#' 
#' 
#' @export
#'
predict_inactivation <- function(times,
                                 primary_model,
                                 environment = "constant",
                                 secondary_models = NULL,
                                 env_conditions = NULL,
                                 ...,
                                 check = TRUE
                                 # logbase_logN = 10,
                                 # formula = . ~ time
) {

  if (environment == "constant") {

    ## Get the name of the primary model and remove it from the vector

    my_model <- primary_model$model
    my_pars <- primary_model
    my_pars$model <- NULL

    ## Undo the log transformations needed in the parameters

    for (each_par in names(my_pars)) {

      if (grepl("log", each_par)) {

        my_pars[[gsub("log", "", each_par)]] <- 10^my_pars[[each_par]]
        my_pars[[each_par]] <- NULL  # To avoid issues later

      }

    }

    ## Give a warning if someone defined environmental conditions

    if (! is.null(env_conditions)) {
      warning("env_conditions are ignored for environment = 'constant'")
    }

    ## Give a warning if someone defined secondary models

    if (! is.null(secondary_models)) {
      warning("secondary_models are ignored for environment = 'constant'")
    }

    ## Call the prediction function

    my_sim <- predict_constant_inactivation(times,
                                  my_model,
                                  my_pars,
                                  check = check
                                  )

    ## Prepare the output class

    out <- list(
      simulation = my_sim,
      primary_model = primary_model,
      pars = my_pars,
      environment = "constant",
      secondary_models = NULL,
      env_conditions = NULL
    )

    ## Calculate the variation of the primary model (constant for isothermal)

    out$primary_pars <- tibble(time = times) %>%
      bind_cols(.,
                as_tibble(my_pars)
      )

    ## Calculate the effects of the secondary models (constant for isothermal)

    sec_effects <- NULL

    out$sec_effects <- sec_effects

    ## Define the class

    class(out) <- c("InactivationPrediction", class(out))

    ## Return

    out

  } else if (environment == "dynamic") {

    ## Check arguments specific for dynamic conditions

    if (check) {

      if (is.null(env_conditions)) {
        stop("env_conditions must be defined for predictions in dynamic environments")
      }

      if (is.null(secondary_models)) {
        stop("secondary_models must be defined for predictions in dynamic environments")
      }

    }

    ## Check that times starts at 0. Give a warning otherwise

    if (min(times) != 0) {
      warning(paste("times does not start at t=0.",
                    "Be mindful that the calculation assumes that the first value of times indicates the initial time for the simulation (i.e., the time point where N = N0)",
                    "If this is not what you intend, just pass c(0, times)."
      )
      )
    }

    ## Make the log transformations needed for the primary model

    my_model <- primary_model$model
    my_pars <- primary_model
    my_pars$model <- NULL

    for (each_par in names(my_pars)) {

      if (grepl("log", each_par)) {

        my_pars[[gsub("log", "", each_par)]] <- 10^my_pars[[each_par]]
        my_pars <- my_pars[names(my_pars) != each_par]  # Remove the log-transformed to avoid issues

      }

    }

    primary_model <- my_pars
    primary_model$model <- my_model

    ## Make the log transformations needed for the secondary model

    for (i in 1:length(secondary_models)) {

      pars <- secondary_models[[i]]

      ## The reference value

      if (grepl("log", pars$par)) {
        secondary_models[[i]]$ref <- 10^pars$ref
        secondary_models[[i]]$par <- gsub("log", "", pars$par)
      }

      ## Parameters of the secondary model for each factor

      pars$par <- NULL
      pars$model <- NULL
      pars$ref <- NULL

      for (each_factor in names(pars)) {

        for (j in 1:length(pars[[each_factor]])) {

          this_par <- names(pars[[each_factor]])[[j]]

          if (grepl("log", this_par)) {

            secondary_models[[i]][[each_factor]][[j]] <- 10^pars[[each_factor]][[j]]
            names(secondary_models[[i]][[each_factor]])[[j]] <- gsub("log", "", this_par)

          }

        }


      }


    }

    ## Call the calculation function

    my_sim <- predict_dynamic_inactivation(times,
                                           primary_model,
                                           env_conditions,
                                           secondary_models,
                                           ...,
                                           check = check
                                           # formula = formula,
                                           )

    ## Prepare the output class

    out <- list(
      simulation = my_sim,
      primary_model = primary_model,
      pars = my_pars,
      environment = "dynamic",
      secondary_models = secondary_models,
      env_conditions = env_conditions
      # logbase_logN = 10,
    )

    ## Calculate the variation of the primary model

    env_interpolator <- approx_env(env_conditions)
    aa <-  map_dfc(env_interpolator, ~.(my_sim$time))

    out$primary_pars <- c(1:nrow(aa)) %>%
      map(
        ~ apply_secondary_models(as.list(aa[.,]), secondary_models)
      ) %>%
      map(
        ~ set_names(.,
                    secondary_models %>% map(~.$par) %>% unlist()
                    )
      ) %>%
      map_dfr(as_tibble) %>%
      mutate(time = my_sim$time) %>%
      select(time, everything())

    ## Calculate the effects on each parameter

    sec_effects <- lapply(secondary_models, function(this_par) {

      c(1:nrow(aa)) %>%
        map_dfr(.,
            ~ get_effects(this_par, aa[.,])
            ) %>%
        mutate(time = my_sim$time) %>%
        select(time, everything())

    })

    names(sec_effects) <- secondary_models %>% map(~.$par) %>% unlist()
    out$sec_effects <- sec_effects


    ## Define the class

    class(out) <- c("InactivationPrediction", class(out))

    ## Return

    out


  } else {

    stop("environment must be 'constant' or 'dynamic', got: ", environment)

  }

}
