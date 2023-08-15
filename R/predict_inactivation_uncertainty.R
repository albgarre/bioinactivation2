
#' AA
#'
#' @importFrom MASS mvrnorm
#' @importFrom rlang .data
#'
#' @export
#'
predict_inactivation_uncertainty <- function(model_name,
                                             times,
                                             n_sims,
                                             pars,
                                             corr_matrix = diag(nrow(pars))
                                             # check = TRUE
) {

  # ## Checks
  #
  # if (isTRUE(check)) {
  #
  #   check_stochastic_pars(model_name, pars, corr_matrix)
  #
  # }

  ## Generate the parameter sample

  mus <- pars$mean
  stdevs <- pars$sd
  b <- stdevs %*% t(stdevs)
  cov_matrix <- b * corr_matrix

  par_sample <- as.data.frame(mvrnorm(n_sims,
                                      mus,
                                      cov_matrix
                                      )
                              ) %>%
    set_names(pars$par)

  ## Undo the transformation

  for (i in 1:nrow(pars)) {
    transf <- pars$scale[i]

    new_col <- switch(transf,
                      original = par_sample[,i],
                      sqrt = par_sample[,i]^2,
                      log = 10^par_sample[,i],
                      stop("Unknown scale:", transf)
    )

    par_sample[,i] <- new_col
  }

  ## Do the simulations

  aa <- par_sample %>%
    mutate(iter = row_number(),
           model = model_name)

  my_sims <- split(aa, aa$iter) %>%
    # split(.$iter) %>%
    map(as.list) %>%
    map(
      ~ predict_inactivation(times, .)
    ) %>%
    imap_dfr(~ mutate(.x$simulation, iter = .y))

  ## Extract the quantiles

  q_values <- my_sims %>%
    group_by(.data$time) %>%
    summarize(q50 = quantile(.data$logN, probs = .5, na.rm=TRUE),
              q10 = quantile(.data$logN, probs = .1, na.rm=TRUE),
              q90 = quantile(.data$logN, probs = .9, na.rm=TRUE),
              q05 = quantile(.data$logN, probs = .05, na.rm=TRUE),
              q95 = quantile(.data$logN, probs = .95, na.rm=TRUE),
              m_logN= mean(.data$logN, na.rm=TRUE)
    )

  ## Prepare the output

  out <- list(
    sample = par_sample,
    simulations = my_sims,
    model = model_name,
    mus = mus,
    sigma = cov_matrix,
    quantiles = q_values
  )

  ## Update the class and return

  class(out) <- c("InactivationUncertainty", class(out))

  out

}


















