
#' Guess for one environmental factor
#'
#' @param model_name a 1D character vector defining the secondary model
#' @param y a numeric vector with the output variable (the parameter values)
#' @param x a numeric vector with the input variable (the environmental factor)
#'
make_guess_factor <- function(model_name, x, y) {

  logy <- log10(y)

  ## Guess for Tref

  xref <- mean(x, na.rm = TRUE)

  ## Guess for z-value

  z <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/(max(logy, na.rm = TRUE) - min(logy, na.rm = TRUE))

  ## Guess for n

  n <- 1

  ## Guess for b-value

  b <- (max(y, na.rm = TRUE) - min(y, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

  ## Guess for Tcrit

  threshold <- min(y, na.rm = TRUE) + .2

  xcrit <- x[y>3][1]

  ## k

  k <- (max(y, na.rm = TRUE) - threshold)/(max(x, na.rm = TRUE) - xcrit)

  ## Guess for Ea

  lny <- log(y)
  invx <- 1/x
  R <- 8.31

  Ea <- R*abs(( max(lny, na.rm = TRUE) - min(lny, na.rm = TRUE) )/( max(invx, na.rm = TRUE) - min(invx, na.rm = TRUE) ))

  ## Select the right parameters

  out <- list(xref = xref, z = z, n = n, xcrit = xcrit, k = k, Ea = Ea, b = b)

  par_map <- tribble(
    ~ par,  ~Bigelow, ~genBigelow, ~Arrhenius, ~Lineal, ~logExponential,
    "xref", TRUE, TRUE, TRUE, TRUE, FALSE,
    "z", TRUE, TRUE, FALSE, FALSE, FALSE,
    "n", FALSE, TRUE, FALSE, FALSE, FALSE,
    "b", FALSE, FALSE, FALSE, TRUE, TRUE,
    "xcrit", FALSE, FALSE, FALSE, FALSE, TRUE,
    "k", FALSE, FALSE, FALSE, FALSE, TRUE,
    "Ea", FALSE, FALSE, TRUE, FALSE, FALSE
  )

  my_pars <- par_map$par[par_map[[model_name]]]
  out <- out[my_pars]

  unlist(out)

}


#' Initial guesses for fitting secondary inactivation models
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function uses some heuristics to provide initial guesses for the parameters
#' of the growth model selected that can be used with [fit_inactivation_secondary()].
#'
#' @export
#'
make_guess_secondary <- function(fit_data,
                                 model_name,
                                 formula = my_par ~ temp
                                 ) {

  # ## Check that we know the model
  #
  # if ( ! (primary_model %in% primary_model_data()) ) {
  #   stop("Unkonwn model: ", primary_model)
  # }
  #

  ## Apply the formula

  y_col <- lhs(formula)

  vars <- all.vars(formula)
  vars <- vars[vars != y_col]

  ## Make the guesses for each factor

  y <- fit_data[[y_col]]

  out <- lapply(vars, function(this_var) {

    ## Get the guess for each factor

    x <- fit_data[[this_var]]

    guess <- make_guess_factor(model_name, x, y)

    names(guess) <- paste0(this_var, "_", names(guess))
    guess

  }) %>%
    unlist()

  ## Add the value at the reference conditions

  c(out, ref = median(y, na.rm = TRUE))

}



#' Initial guesses for fitting primary inactivation models
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function uses some heuristics to provide initial guesses for the parameters
#' of the growth model selected that can be used with [fit_inactivation()].
#'
#' @param fit_data the experimental data. A tibble (or data.frame) with a column
#' named `time` with the elapsed time and one called `logN` with the logarithm
#' of the population size
#' @param primary_model a string defining the equation of the primary model,
#' as defined in [primary_model_data()]
#' @param formula an object of class "formula" describing the x and y variables.
#' `logN ~ time` as a default.
#'
#' @return A named numeric vector of initial guesses for the model parameters
#'
#' @importFrom tibble tribble
#'
#' @export
#'
#'
make_guess_primary <- function(fit_data, primary_model,
                               formula = logN ~ time
) {

  ## Check that we know the model

  if ( ! (primary_model %in% primary_model_data()) ) {
    stop("Unkonwn model: ", primary_model)
  }

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

  if (primary_model == "Weibull_2phase") {
    stop("Automagicical initial guesses cannot be calculated for the 2-phase Weibull model")
  }

  ## Guess for logN0

  logN0 <- max(fit_data$logN, na.rm = TRUE)

  ## Guess for SL

  SL <- min(fit_data$time[which(fit_data$logN < logN0 - .5)], na.rm = TRUE)

  ## Gues for logNres

  logNres <- min(fit_data$logN, na.rm = TRUE)

  ## Guess for logD

  tmax <- max(fit_data$time, na.rm = TRUE)
  D <- tmax/(logN0 - logNres)
  logD <- log10(D)

  ## Guess for logdelta

  logdelta <- logD

  ## Guess for p

  p <- 1


  ## Guess for n

  n <- 1

  ## Guess for b

  b <- 1/D

  ## Guess for k

  k <- 1/D

  ## Guess for Delta

  Delta <- 5

  ## return

  out <- list(logN0 = logN0, SL = SL, logNres = logNres,
              logD = logD, p = p, n = n, b = b, k = k,
              logdelta = logdelta)

  par_map <- tribble(
    ~ par,  ~Bigelow, ~Peleg, ~Mafart, ~Geeraerd, ~Metselaar, ~Trilinear, ~Geeraerd_k, ~Geeraerd_k_noTail, ~Geeraerd_noTail, ~Geeraerd_k_noShoulder, ~Geeraerd_noShoulder,
    "logN0", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    "SL", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
    "logNres", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE,
    "logD", TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE,
    "p", FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    "logdelta", FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    "n", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    "b", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    "k", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE,
    "Delta", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
  )

  my_pars <- par_map$par[par_map[[primary_model]]]

  unlist(out[my_pars])

}









