
#'
#'
get_effects <- function(this_par, env_conditions) {

  ## Extract the info

  par_name <- this_par$par
  this_model <- this_par$model
  ref_value <- this_par$ref

  this_par$par <- NULL
  this_par$model <- NULL
  this_par$ref <- NULL

  ## Calculate the effect for each factor

  # browser()

  effects <- this_par %>%
    imap_dbl(.,
             ~ switch(this_model,
                      Bigelow = sec_Bigelow(
                        x = env_conditions[[.y]],
                        xref = .x[["xref"]], z = .x[["z"]]
                      ),
                      genBigelow = sec_genBigelow(
                        x = env_conditions[[.y]],
                        xref = .x[["xref"]], z = .x[["z"]], n = .x[["n"]]
                      ),
                      Lineal = sec_Lineal(
                        x = env_conditions[[.y]],
                        xref = .x[["xref"]], b = .x[["b"]]
                      ),
                      Arhenius = sec_Arhenius(
                        x = env_conditions[[.y]],
                        xref = .x[["xref"]], Ea = .x[["Ea"]]
                      ),
                      LogExponential = sec_logExponential(
                        x = env_conditions[[.y]],
                        k = .x[["k"]], xc = .x[["Xcrit"]]),
                      stop(paste("Unknown model:", this_model))
             )
             # ~ sec_Bigelow(x = env_conditions[[.y]], xref = .x[["xref"]], z = .x[["z"]])
    )

  effects

}


#' Calculates the effect of the secondary models on the model parameters for some environmental conditions
#'
#' @param env_conditions A list with as many entries as environmental factors. Each is
#' a numeric vector of length 1 (a tibble of 1 row also works).
#' @param secondary_models A nested list describing the secondary models (as in [predict_inactivation()]).
#'
#' @returns A list, with as many entries as model parameters, with the parameter values.
#'
apply_secondary_models <- function(env_conditions, secondary_models) {

  out <- lapply(secondary_models, function(this_par) {

    effects <- get_effects(this_par, env_conditions)

    ## Put together and apply

    full_effect <- sum(effects)

    if (this_par$model %in% c("Bigelow", "genBigelow")) {

      10^(log10(this_par$ref) + full_effect)

    } else {

      this_par$ref + full_effect

    }

  })

  out

}

#' Generates functions for linear interpolation of environmental conditions
#'
#' @param env_conditions A tibble describing the variation of the environmental
#' conditions through the storage time. Must contain a column named `time`
#' and as many additional columns as environmental factors.
#'
#' @importFrom stats approxfun
#'
#' @return A list of functions that return the value of each environmental
#' condition for some storage time
#'
approx_env <- function(env_conditions) {

  out <- lapply(names(env_conditions[-1]), function(this_col) {

    x <- env_conditions$time
    y <- env_conditions[[this_col]]

    approxfun(x, y, rule = 2)

  })

  names(out) <- names(env_conditions[-1])
  out

}



#' A helper to convert guesses to the format for predictions
#'
#' @importFrom stringr str_detect str_replace
#'
#'
convert_dynamic_guess <- function(sec_models, guess, known) {

  p <- c(guess, known)
  p_names <- names(p)

  out <- list()

  for (each_model in sec_models) {

    ## Type of secondary value

    this_model <- list(
      model = each_model$model
    )

    ## Reference value

    this_par <- each_model$par

    if (paste0(this_par, "_logref") %in% p_names) {

      this_model$par <- paste0("log", this_par)
      this_model$ref <- p[[paste0(this_par, "_logref")]]

    } else {
      this_model$par <- this_par
      this_model$ref <- p[[paste0(this_par, "_ref")]]
    }

    ## Parameters for each factor

    # browser()

    for (each_factor in each_model$depends_on) {

      # browser()

      start <- paste0(this_par, "_", each_factor, "_")

      this_pars <- p[str_detect(names(p), start)]
      names(this_pars) <- str_replace(names(this_pars), start, "")
      this_model[[each_factor]] <- this_pars

    }

    ## Append the model

    out[[each_model$par]] <- this_model

  }

  ## Return

  out

}

