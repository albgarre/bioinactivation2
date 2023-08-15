
#'
#' @export
#'
predict_inactivation <- function(times,
                                 primary_model,
                                 environment = "constant",
                                 secondary_models = NULL,
                                 env_conditions = NULL
                                 # ...,
                                 # check = TRUE,
                                 # logbase_mu = logbase_logN,
                                 # logbase_logN = 10,
                                 # formula = . ~ time
) {

  if (environment == "constant") {

    ## Get the name of the primary model and remove it from the vector

    my_model <- primary_model$model
    my_pars <- primary_model
    my_pars$model <- NULL

    ## Make the log transformations needed

    for (each_par in names(my_pars)) {

      if (grepl("log", each_par)) {

        my_pars[[gsub("log", "", each_par)]] <- 10^my_pars[[each_par]]

      }

    }

    # ## Make it work both with N0 and logN0
    #
    # if ("N0" %in% names(my_pars)) {
    #   my_pars$logN0 <- log10(my_pars$N0)
    # }

    # if (check) {
    #
    #   if (is.null(my_model)) {
    #     stop("primary model must include a 'model' entry")
    #   }
    #
    # }
#
#     ## Give a warning if someone defined environmental conditions
#
#     if (! is.null(env_conditions)) {
#       warning("env_conditions are ignored for 'constant' predictions")
#     }
#
#     ## Give a warning if someone defined secondary models
#
#     if (! is.null(secondary_models)) {
#       warning("secondary_models are ignored for 'constant' predictions")
#     }

    ## Call the prediction function

    my_sim <- predict_constant_inactivation(times,
                                  my_model,
                                  my_pars #,
                                  # ...,
                                  # check = TRUE,
                                  # logbase_logN = 10
                                  )

    ## Prepare the output class

    out <- list(
      simulation = my_sim,
      primary_model = primary_model,
      pars = my_pars,
      environment = "constant",
      secondary_models = NULL,
      env_conditions = NULL
      # logbase_logN = 10,
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

    # if (check) {
    #
    #   if (is.null(env_conditions)) {
    #     stop("env_conditions must be defined for predictions in dynamic environments")
    #   }
    #
    #   if (is.null(secondary_models)) {
    #     stop("secondary_models must be defined for predictions in dynamic environments")
    #   }
    #
    # }
    #
    # ## Check that times starts at 0. Give a warning otherwise
    #
    # if (min(times) != 0) {
    #   warning(paste("times does not start at t=0.",
    #                 "Be mindful that the calculation assumes that the first value of times indicates the initial time for the simulation (i.e., the time point where N = N0)",
    #                 "If this is not what you intend, just pass c(0, times)."
    #   )
    #   )
    # }

    ## Make the log transformations needed for the primary model

    my_model <- primary_model$model
    my_pars <- primary_model
    my_pars$model <- NULL

    for (each_par in names(my_pars)) {

      if (grepl("log", each_par)) {

        my_pars[[gsub("log", "", each_par)]] <- 10^my_pars[[each_par]]

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
                                           secondary_models
                                           # ... #,
                                           # check = check,
                                           # formula = formula,
                                           # logbase_mu = logbase_mu,
                                           # logbase_logN = logbase_logN
                                           )

    # ## Return
    #
    # out

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
