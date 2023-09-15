
#' InactivationFit class
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' The `InactivationFit` class contains an inactivation model fitted to data under
#' static or dynamic conditions. Its constructor is [fit_inactivation()].
#'
#' It is a subclass of list with the items:
#'
#' - approach: fitting approach as in [fit_inactivation()]
#' - algorithm: type of algorithm as in [fit_inactivation()]
#' - data: data used for model fitting
#' - guess: initial guess of the model parameters
#' - known: fixed model parameters
#' - primary_model: a character describing the primary model
#' - fit_results: an instance of modFit or modMCMC with the results of the fit
#' - best_prediction: Instance of [InactivationPrediction] with the best inactivation fit
#' - sec_models: a named vector with the secondary models assigned for each
#' environmental factor. `NULL` for `environment="constant"`
#' - env_conditions: a tibble with the environmental conditions used for model
#' fitting. `NULL` for `environment="constant"`
#' - niter: number of iterations of the Markov chain. `NULL` if `algorithm != "MCMC"`
#' - logbase_logN: base of the logarithm for the definition of the population size
#' (check the relevant vignette)
#' - approach_logN0: Approach for the initial microbial count when fitting several experiments
#'
#' @name InactivationFit
#'
NULL

#' @describeIn InactivationFit vector of fitted model parameters.
#'
#' @param object an instance of [InactivationFit].
#' @param ... ignored
#'
#' @importFrom stats coef
#'
#' @export
#'
coef.InactivationFit <- function(object, ..., step = 2) {

  if (object$approach == "two-steps") {

    if (step == 2) {

      my_models <- object$fit_results$secondary

      out <- lapply(my_models, function(this_model) {

        pars <- coef(this_model)

        if (!is.null(pars)) {
          names(pars) <- paste0(this_model$par_name, "_", names(pars))
        }

        pars

      })

      names(out) <- NULL
      unlist(out)

    } else if (step == 1) {

      object$primary_pars

    } else {
      stop(paste("'step' must be 1 or 2, got:", step))
    }

  } else if (object$algorithm == "regression") {

    coef(object$fit_results)

  } else {
    object$fit_results$bestpar
  }

}

#' @describeIn InactivationFit statistical summary of the fit.
#'
#' @param object Instance of [InactivationFit]
#' @param ... ignored
#'
#' @export
#'
summary.InactivationFit <- function(object, ..., step = 2) {

  if (object$approach == "two-steps") {

    if (step == 2) {

      out <- object$fit_results$secondary %>% map(summary)

    } else if (step == 1) {

      out <- object$fit_results$primary %>% map(summary)

    } else {

      stop(paste("'step' must be 1 or 2, got:", step))

    }
  } else {

    out <- summary(object$fit_results)

  }

  out

}

#' @describeIn InactivationFit vector of model predictions.
#'
#' @param object an instance of [InactivationFit]
#' @param ... ignored
#' @param times numeric vector describing the time points for the prediction.
#' If `NULL` (default), uses the same points as those used for fitting.
#' @param env_conditions tibble describing the environmental conditions as in [fit_inactivation()].
#' If `NULL` (default), uses the environmental condition of the fitting.
#'
#' @export
#'
predict.InactivationFit <- function(object, times = NULL, env_conditions = NULL, ...) {
  
  if (is.null(times)) {  ## Used the times of the data if NULL

    times <- object$data$time

  }

  if (object$approach == "primary") {  # Prediction under constant environment

    pars <- c(coef(object), unlist(object$known))
    my_model <- as.list(pars)
    my_model$model <- object$primary_model

    pred <- predict_inactivation(times, my_model
                                 #, check = FALSE,
                                 # logbase_logN = object$logbase_logN
                                 )
    pred$simulation$logN

  } else if (object$approach == "dynamic") {  ## Prediction under dynamic conditions

    if (is.null(env_conditions)) {  # Used the environment of the data if NULL

      env_conditions <- object$env_conditions

    }

    pred <- predict_inactivation(environment = "dynamic",
                                 times,
                                 object$best_prediction$primary_model,
                                 object$best_prediction$secondary_models,
                                 env_conditions
                                 )
    
    pred$simulation$logN

  } else if (object$approach == "one-step") {  # Predictions using a one-step model

    ## Make the primary model

    # p <- c(coef(object), object$known)
    #
    # primary_model <- list(model = object$primary_model)
    #
    # initial <- p[str_detect(names(p), "N0")]
    # aa <- initial
    # names(aa) <- NULL
    # primary_model[[names(initial)]] <- aa
    #
    # if (str_detect(primary_model_name, "Geeraerd")) {
    #   initial <- p[str_detect(names(p), "C0")]
    #   aa <- initial
    #   names(aa) <- NULL
    #   primary_model[[names(initial)]] <- aa
    # }
    #
    #
    # pred <- predict_inactivation(environment = "dynamic",
    #                              times,
    #                              object$best_prediction$primary_model,
    #                              object$sec_models,
    #                              env_conditions
    #                              # logbase_logN = object$logbase_logN
    # )
    #
    # pred$simulation$logN
    stop("predict method not implemented for this approach")

  } else if (object$approach == "global") {
    
    if ( is.null(times) & is.null(env_conditions) ) {
      fitted(object)
    } else {
      
      pred <- predict_inactivation(environment = "dynamic",
                                   times,
                                   object$best_prediction[[1]]$primary_model,
                                   object$best_prediction[[1]]$secondary_models,
                                   env_conditions
                                   )
      
      pred$simulation$logN
      
    }
    
  } else {
    stop("predict method not implemented for this approach")
  }


}

#' @describeIn InactivationFit vector of model residuals.
#'
#' @param object Instance of [InactivationFit]
#' @param ... ignored
#'
#' @importFrom stats residuals
#'
#' @export
#'
residuals.InactivationFit <- function(object, ...) {

  if (object$algorithm == "MCMC") {

    pred <- predict(object)

    pred - object$data$logN

  } else {

    residuals(object$fit_results)

  }

}

#' @describeIn InactivationFit variance-covariance matrix of the model, estimated
#' as 1/(0.5*Hessian) for regression and as the variance-covariance of the draws
#' for MCMC
#'
#' @param object an instance of [InactivationFit]
#' @param ... ignored
#'
#' @export
#'
vcov.InactivationFit <- function(object, ...) {

  if (object$algorithm == "MCMC") {

    cov(object$fit_results$pars)

  } else {

    # The code has been adapted from the one of summary.modFit

    covar  <- try(solve(0.5*object$fit_results$hessian), silent = TRUE)

    if (!is.numeric(covar)) {
      warning("Cannot estimate covariance; system is singular")

      param  <- object$par
      p      <- length(param)

      covar <- matrix(data = NA, nrow = p, ncol = p)
    }

    covar

  }
}

#' @describeIn InactivationFit deviance of the model.
#'
#' @param object an instance of [InactivationFit]
#' @param ... ignored
#'
#' @importFrom stats deviance
#'
#' @export
#'
deviance.InactivationFit <- function(object, ...) {

  if (object$algorithm == "MCMC") {

    sum(residuals(object)^2)

  } else {
    deviance(object$fit_results)
  }

}

#' @describeIn InactivationFit vector of fitted values.
#'
#' @param object an instance of [InactivationFit]
#' @param ... ignored
#' 
#' @importFrom purrr imap_dfr
#' @importFrom tidyselect everything
#' 
#'
#' @export
#'
fitted.InactivationFit <- function(object, ...) {
  
  if (object$approach == "global") {
    
    object$data %>%
      imap_dfr(~ mutate(.x, exp = .y)) %>%
      select(exp, everything()) %>%
      mutate(res = residuals(object),
             fitted = res + logN) %>%
      pull(fitted)

  } else if (object$approach == "one-step") {

    ## Make the primary model

    p <- c(coef(object), object$known)

    primary_model <- list(model = object$primary_model)

    initial <- p[str_detect(names(p), "N0")]
    aa <- initial
    names(aa) <- NULL
    primary_model[[names(initial)]] <- aa

    if (str_detect(object$primary_model, "Geeraerd")) {
      initial <- p[str_detect(names(p), "C0")]
      aa <- initial
      names(aa) <- NULL
      primary_model[[names(initial)]] <- aa
    }

    ## Make the predictions

    cond <- select(object$data, -logN)

    pred <- lapply(1:nrow(cond), function(i) {

      env <- bind_rows(cond[i,], cond[i,])
      env$time[1] <- 0

      predict_inactivation(environment = "dynamic",
                           seq(0, env$time[2], length = 10),
                           primary_model,
                           object$sec_models,
                           env)$simulation %>%
        tail(1) %>%
        pull(logN)

    }) %>%
      unlist()

    pred

  } else {

    predict(object)

  }
}

#' @describeIn InactivationFit loglikelihood of the model
#'
#' @param object an instance of InactivationFit
#' @param ... ignored
#'
#' @export
#'
logLik.InactivationFit <- function(object, ...) {
  
  ## Get the number of data points
  
  if (object$approach == "global") {
    n <- object$data %>% map_dbl(nrow) %>% sum()
  } else {
    n <- nrow(object$data)
  }

  if (object$algorithm == "regression") {

    # n <- nrow(object$data)
    sigma <- sqrt(object$fit_results$ssr/object$fit_results$df.residual)

    lL <- - n/2*log(2*pi) -n/2 * log(sigma^2) - 1/2/sigma^2*object$fit_results$ssr

    lL

  } else {

    # n <- nrow(object$data)
    SS <- min(object$fit_results$SS, na.rm = TRUE)

    df <- n - length(coef(object))

    sigma <- sqrt(SS/df)

    lL <- - n/2*log(2*pi) -n/2 * log(sigma^2) - 1/2/sigma^2*SS

    lL

  }
}

#' @describeIn InactivationFit Akaike Information Criterion
#'
#' @param object an instance of InactivationFit
#' @param ... ignored
#' @param k penalty for the parameters (k=2 by default)
#'
#' @importFrom stats logLik
#'
#' @export
#'
AIC.InactivationFit <- function(object, ..., k=2) {
  
  ## Get the number of data points
  
  if (object$approach == "global") {
    n <- object$data %>% map_dbl(nrow) %>% sum()
  } else {
    n <- nrow(object$data)
  }

  ## Normal AIC

  p <- length(coef(object))

  lL <- logLik(object)

  AIC <- 2*p - 2*lL

  ## Calculate the penalty

  # n <- nrow(object$data)

  penalty <- (k*p^2 + k*p)/(n - p - 1)

  ## Return

  AIC + penalty

}

#' @describeIn InactivationFit compares the fitted model against the data.
#'
#' @param x The object of class [InactivationFit] to plot.
#' @param y ignored
#' @param ... ignored.
#' @param add_factor whether to plot also one environmental factor.
#' If `NULL` (default), no environmental factor is plotted. If set
#' to one character string that matches one entry of x$env_conditions,
#' that condition is plotted in the secondary axis. Ignored if `environment="constant"`
#' @param ylims A two dimensional vector with the limits of the primary y-axis.
#' `NULL` by default
#' @param label_y1 Label of the primary y-axis.
#' @param label_y2 Label of the secondary y-axis. Ignored if `environment="constant"`
#' @param line_col Aesthetic parameter to change the colour of the line geom in the plot, see: [geom_line()]
#' @param line_size Aesthetic parameter to change the thickness of the line geom in the plot, see: [geom_line()]
#' @param line_type Aesthetic parameter to change the type of the line geom in the plot, takes numbers (1-6) or strings ("solid") see: [geom_line()]
#' @param point_col Aesthetic parameter to change the colour of the point geom, see: [geom_point()]
#' @param point_size Aesthetic parameter to change the size of the point geom, see: [geom_point()]
#' @param point_shape Aesthetic parameter to change the shape of the point geom, see: [geom_point()]
#' @param line_col2 Same as lin_col, but for the environmental factor. Ignored if `environment="constant"`
#' @param line_size2 Same as line_size, but for the environmental factor. Ignored if `environment="constant"`
#' @param line_type2 Same as lin_type, but for the environmental factor. Ignored if `environment="constant"`
#' @param label_x Label of the x-axis
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_point
#' @importFrom rlang .data
#' @importFrom graphics plot
#' @importFrom cowplot theme_cowplot plot_grid
#'
plot.InactivationFit <- function(x, y=NULL, ...,
                           add_factor = NULL,
                           line_col = "black",
                           line_size = 1,
                           line_type = 1,
                           point_col = "black",
                           point_size = 3,
                           point_shape = 16,
                           ylims = NULL,
                           label_y1 = NULL,
                           label_y2 = add_factor,
                           label_x = "time",
                           line_col2 = "black",
                           line_size2 = 1,
                           line_type2 = "dashed",
                           type = 1) {

  ## Get the label for the y-axis

  # logbase <- x$logbase_logN
  #
  # if ( abs(logbase - exp(1)) < .1 ) {
  #   logbase <- "e"
  # }

  if (is.null(label_y1)) {
    # label_y1 <- paste0("logN (in log-", logbase, ")")
    label_y1 <- "logN"
  } else {
    label_y1 <- label_y1
  }

  if (x$approach %in% c("primary")) {

    p <- plot(x$best_prediction,
              line_col = line_col,
              line_size = line_size,
              line_type = line_type,
              ylims = ylims,
              label_y1 = label_y1,
              label_x = label_x)

    p + geom_point(aes(x = .data$time, y = .data$logN), data = x$data,
                   col = point_col,  size = point_size, shape = point_shape) +
      theme_cowplot()

  } else if (x$approach == "dynamic") {

    p <- plot(x$best_prediction
              # add_factor = add_factor,
              # ylims = ylims,
              # label_y1 = label_y1,
              # label_y2 = label_y2,
              # line_col = line_col,
              # line_size = line_size,
              # line_type = line_type,
              # line_col2 = line_col2,
              # line_size2 = line_size2,
              # line_type2 = line_type2,
              # label_x = label_x
    )

    p + geom_point(aes(x = .data$time, y = .data$logN), data = x$data,
                   col = point_col,  size = point_size, shape = point_shape) +
      theme_cowplot()

  } else if(x$approach == "one-step") {

    if (type == 1) {  # predicted versus observed

      x$data %>%
        mutate(pred = fitted(x)) %>%
        ggplot(aes(x = logN, y = pred)) +
        geom_point() +
        geom_smooth(approach = "lm", se = FALSE) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "gray")

    } else if (type == 2) {  # facet-based

      ## TODO

    } else {

      stop(paste0("type must be 1 or 2 for approach = 'one-step'. Got: ", type))

    }
  } else if (x$approach == "global") {
    
    p <- x$best_prediction %>% map(plot)
    
    p <- map2(p, x$data,
         ~ .x + geom_point(aes(x = time, y = logN), data = .y)
         )
    
    plot_grid(plotlist = p, labels = names(p))
    
  } else if (x$approach == "two-steps") {
    
    if (type == 1) {
      
      p <- x$fit_results$primary %>%
        map( ~ plot(.) )
      
      plot_grid(plotlist = p, labels = names(p))

    } else {
      stop(paste("type must be 1 for this approach, got:", type))
    }
    
    
  } else {
    stop("plot() method not implemented for this approach")
  }

}




























