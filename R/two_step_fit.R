
#' AA
#'
#'
fit_two_step <- function(fit_data,
                         model_name,
                         start,
                         known,
                         upper,
                         lower,
                         secondary_models
                         # algorithm,
                         # niter
                         # approach_logN0 = approach_logN0
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

  ## Split the data

  split_data <- fit_data %>%
    unite(condition, -c(time, logN), sep = "_", remove = FALSE) %>%
    split(.$condition)

  ## Initial guesses for the primary fit

  initial_guesses <- split_data %>%
    map(
      ~ make_guess_primary(., model_name)
    )

  # ## Remove the known parameters from the guesses
  #
  # if (length(known) > 0) {
  #
  #   initial_guesses <- as.list(initial_guesses)
  #
  #   known_pars <- names(known)
  #
  #   for (each_par in known_pars) {
  #
  #     initial_guesses <- initial_guesses %>%
  #       map(~ remove_stuff(., each_par))
  #   }
  #
  #   unlist(initial_guesses)
  #
  # }

  ## Fit the primary models

  primary_fits <- map2(initial_guesses, split_data,
                       ~ fit_inactivation("primary",
                                          select(.y, time, logN),
                                          model_name,
                                          .x, 
                                          known,
                                          check = FALSE
                                          )
                       )

  ## Extract the parameters

  my_pars <- primary_fits %>% map(coef) %>% map(as.list)

  ## Join them with the data

  fitted_pars <- split_data %>%
    map( ~ select(., -time, -logN)) %>%
    map( ~ head(., 1)) %>%
    map2_dfr(my_pars, bind_cols) %>%
    select(-condition)

  ## Undo the logtransformations

  par_names <- names(fitted_pars)

  for (each_name in par_names) {

    if (grepl("log", each_name)) {

      fitted_pars[[str_replace(each_name, "log", "")]] <- 10^fitted_pars[[each_name]]
      fitted_pars <- select(fitted_pars, - matches(each_name))

    }

  }

  ## Do the 2nd step fit

  secondary_fits <- lapply(secondary_models, function(each_par) {

    if (length(each_par$depends_on) == 0) {  # No secondary model to fit

      NULL

    } else {

      ## Extract the data we are using

      this_d <- fitted_pars %>%
        select(matches(each_par$depends_on), matches(each_par$par))

      ## Get initial guesses

      my_formula <- as.formula(
        paste(each_par$par, " ~ ",
              paste(each_par$depends_on, collapse = " + ")
              )
      )

      my_guess <- make_guess_secondary(this_d, each_par$model,
                                       formula = my_formula
                                       )

      # ## Remove the known parameters from the guesses

      this_known <- known[grepl(each_par$par, names(known))]

      if (length(known) > 0) {

        ## Remove the part of the name we dont need

        names(this_known) <- str_replace(names(this_known), paste0(each_par$par, "_"), "")

        my_guess <- my_guess[ ! (names(my_guess) %in% names(this_known)) ]

      }

      # browser()

      ## Fit the secondary model

      fit_inactivation_secondary(this_d,
                                 each_par$model,
                                 guess = my_guess,
                                 known = this_known,
                                 formula = my_formula
                                 )

    }

  })

  names(secondary_fits) <- secondary_models %>% map(~ .$par)

  ## Prepare the output

  out <- list(
    fit_data = fit_data,
    model_name = model_name,
    start = start,
    known = known,
    # upper = upper,
    # lower = lower,
    sec_models = secondary_models,
    fit_results = list(
      primary = primary_fits,
      secondary = secondary_fits
    ),
    primary_pars = fitted_pars,
    algorithm = "regression"

  )

  ## Return

  out

}

#'
#'
remove_stuff <- function(x, stuff) {

  x[[stuff]] <- NULL

}



