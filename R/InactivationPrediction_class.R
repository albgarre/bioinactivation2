
#' InactivationPrediction class
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' The `InactivationPrediction` class contains the results of an inactivation prediction.
#' Its constructor is [predict_inactivation()].
#'
#' It is a subclass of list with the items:
#'
#' - simulation: a tibble with the model simulation
#' - primary model: a list describing the primary model as in [predict_inactivation()]
#' - environment: a character describing the type of environmental conditions
#' as in [predict_inactivation()]
#' - env_conditions: a named list with the functions used to approximate the (dynamic)
#' environmental conditions. `NULL` if `environment="constant"`.
#' - secondary_models: a named list describing the secondary models as in [predict_inactivation()].
#' `NULL` if `environment="constant"`.
#' - sec_effects: a tibble describing the effect of each environmental factor per time point.
#' - logbase_mu: the log-base for the definition of parameter mu (see the relevant vignette)
#' - logbase_logN: the log-base for the definition of the logarithm of the population size
#'
#' @name InactivationPrediction
#'
NULL

#' @describeIn InactivationPrediction print of the model
#'
#' @param x An instance of `InactivationPrediction`.
#' @param ... ignored
#'
#' @export
#'
print.InactivationPrediction <- function(x, ...) {


  if (x$environment == "constant") {

    cat("Inactivation prediction based on primary models\n\n")

    cat(paste("Inactivation model:", x$primary_model$model, "\n\n"))

    cat("Parameters of the primary model:\n")
    print(coef(x))

    # logbase <- x$logbase_logN
    #
    # if ( abs(logbase - exp(1)) < .1 ) {
    #   logbase <- "e"
    # }
    # cat("\n")
    # cat(paste0("Population size defined in log-", logbase, " scale\n"))

  } else if (x$environment == "dynamic") {

    cat("Inactivation prediction under dynamic environmental conditions\n\n")

    env <- names(x$env_conditions)
    cat(paste("Environmental factors included:", paste(env, collapse = ", "), "\n\n"))

    cat("Primary model: ")
    cat(x$primary_model$model)
    cat("\n")

    cat("Parameters of the primary model:\n")
    p <- x$primary_model
    p$model <- NULL
    print(unlist(p))
    cat("\n")

    # cat(paste0("Population size defined in log-", logbase, " scale\n\n"))
    #
    for (i in 1:length(x$secondary_models)) {
      cat(paste("Secondary model for parameter '", x$secondary_models[[i]]$par, "': ",
                x$secondary_models[[i]]$model, "\n", sep = ""))

      p <- x$secondary_models[[i]]
      p$model <- NULL
      p$par <- NULL

      cat(paste("Value under reference conditions: ", p$ref, "\n",
                sep = ""))
      p$ref <- NULL

      for (each_factor in names(p)) {

        cat(paste("Parameters with respect to ", each_factor, ":\n",
                  sep = ""))
        print(p[[each_factor]])
        cat("\n")

      }

      cat("\n")
    }

  }

}

#' @describeIn InactivationPrediction coefficients of the model
#'
#' @param object an instance of [InactivationPrediction]
#' @param ... ignored
#'
#' @export
#'
coef.InactivationPrediction <- function(object, ...) {

  if (object$environment == "constant") {

    out <- object$primary_model
    out$model <- NULL
    unlist(out)

  } else if (object$environment == "dynamic") {

    list(
      primary = object$primary_model,
      secondary = object$secondary_models
    )

  }

}


#' @describeIn InactivationPrediction predicted inactivation curve.
#'
#' @param x The object of class `InactivationPrediction` to plot.
#' @param y ignored
#' @param ... ignored
#' @param add_factor whether to plot also one environmental factor.
#' If `NULL` (default), no environmental factor is plotted. If set
#' to one character string that matches one entry of x$env_conditions,
#' that condition is plotted in the secondary axis. Ignored for `environment="constant"`.
#' @param ylims A two dimensional vector with the limits of the primary y-axis.
#' @param label_y1 Label of the primary y-axis.
#' @param label_y2 Label of the secondary y-axis.
#' @param line_col Aesthetic parameter to change the colour of the line geom in the plot, see: [geom_line()]
#' @param line_size Aesthetic parameter to change the thickness of the line geom in the plot, see: [geom_line()]
#' @param line_type Aesthetic parameter to change the type of the line geom in the plot, takes numbers (1-6) or strings ("solid") see: [geom_line()]
#' @param line_col2 Same as lin_col, but for the environmental factor.
#' @param line_size2 Same as line_size, but for the environmental factor.
#' @param line_type2 Same as lin_type, but for the environmental factor.
#' @param label_x Label of the x-axis.
#' @param type type of plot. See details.
#'
#' @export
#'
plot.InactivationPrediction <- function(x, y=NULL, ...,
                                  add_factor = NULL,
                                  ylims = NULL,
                                  label_y1 = NULL,
                                  label_y2 = add_factor,
                                  line_col = "black",
                                  line_size = 1,
                                  line_type = "solid",
                                  line_col2 = "black",
                                  line_size2 = 1,
                                  line_type2 = "dashed",
                                  label_x = "time",
                                  type = 1
) {

  ## Get the label for the y-axis

  # logbase <- x$logbase_logN
  #
  # if ( abs(logbase - exp(1)) < .1 ) {
  #   logbase <- "e"
  # }
  #
  if (is.null(label_y1)) {
    # label_y1 <- paste0("logN (in log-", logbase, ")")
    label_y1 <- "logN"
  } else {
    label_y1 <- label_y1
  }

  switch(x$environment,
         constant = plot_constant_inactivation(x,
                                          line_col = line_col,
                                          line_size = line_size,
                                          line_type = line_type,
                                          ylims = ylims,
                                          label_y = label_y1,
                                          label_x = label_x
         ),
         dynamic = plot_dynamic_inactivation(x,
                                             add_factor = add_factor,
                                             ylims = ylims,
                                             label_y1 = label_y1,
                                             label_y2 = label_y2,
                                             line_col = line_col,
                                             line_size = line_size,
                                             line_type = line_type,
                                             line_col2 = line_col2,
                                             line_size2 = line_size2,
                                             line_type2 = line_type2,
                                             label_x = label_x,
                                             type = type
         )
  )


}

#' Plotting of inactivation under constant conditions
#' 
#' @importFrom ggplot2 scale_y_continuous xlab 
#' @importFrom rlang .data
#'
plot_constant_inactivation <- function(x, y=NULL, ...,
                                 line_col = "black",
                                 line_size = 1,
                                 line_type = "solid",
                                 ylims = NULL,
                                 label_y = NULL,
                                 label_x = "time") {

  ggplot(x$simulation) +
    geom_line(aes(x = .data$time, y = .data$logN),
              col = line_col,
              size = line_size,
              linetype = line_type) +
    scale_y_continuous(limits = ylims,
                       name = label_y) +
    xlab(label_x) +
    theme_cowplot()

}

#' Plots for fynamic inactivation experiments
#' 
#' @importFrom cowplot plot_grid
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#'
plot_dynamic_inactivation <- function(x,
                                      add_factor = add_factor,
                                      ylims = ylims,
                                      label_y1 = label_y1,
                                      label_y2 = label_y2,
                                      line_col = line_col,
                                      line_size = line_size,
                                      line_type = line_type,
                                      line_col2 = line_col2,
                                      line_size2 = line_size2,
                                      line_type2 = line_type2,
                                      label_x = label_x,
                                      type = type
) {

  if (type == 1) {  # Plot of the survivor curve

    x$simulation %>%
      ggplot() + geom_line(aes(.data$time, .data$logN))

  } else if (type == 2) {  # Plot of the variation of the parameters

    x$primary_pars %>%
      pivot_longer(-"time") %>%
      ggplot() +
      geom_line(aes(x = .data$time, y = .data$value, colour = .data$name)) +
      facet_wrap("name", scales = "free")

  } else if (type == 3) {  # Plot of the effects of each factor
    
    my_plots <- lapply(x$sec_effects, function(this_par) {
      
      if (ncol(this_par) == 1) {  # Constant
        
        ggplot()
        
      } else {
        
        this_par %>%
          pivot_longer(-"time") %>%
          ggplot() +
          geom_line(aes(x = .data$time, y = .data$value, colour = .data$name))
        
      }
      
    })
    
    plot_grid(plotlist = my_plots,
              labels = names(x$sec_effects)
              )

  } else {
    stop(paste0("type must be 1, 2 or 3. Got: ", type))
  }

}


#' @describeIn InactivationPrediction summary of the model
#'
#' @param object An instance of `InactivationPrediction`.
#' @param ... ignored
#'
#' @export
#'
summary.InactivationPrediction <- function(object, ...) {


  if (object$environment == "constant") {

    cat("Inactivation prediction based on primary models\n\n")

    cat(paste("Inactivation model:", object$primary_model$model, "\n\n"))

    cat("Parameters of the primary model:\n")
    print(coef(object))

    # logbase <- object$logbase_logN
    #
    # if ( abs(logbase - exp(1)) < .1 ) {
    #   logbase <- "e"
    # }
    # cat("\n")
    # cat(paste0("Population size defined in log-", logbase, " scale\n"))

    cat("\n")
    cat(paste0("Maximum elapsed time: ", max(object$simulation$time, na.rm=TRUE), "\n"))
    cat(paste0("Population size at the end: ", min(object$simulation$logN, na.rm=TRUE), "\n"))

  } else if (object$environment == "dynamic") {

    cat("Inactivation prediction under dynamic environmental conditions\n\n")

    env <- names(object$env_conditions)
    cat(paste("Environmental factors included:", paste(env, collapse = ", "), "\n\n"))

    cat("Primary model: ")
    cat(object$primary_model$model)
    cat("\n")

    cat("Parameters of the primary model:\n")
    p <- object$primary_model
    p$model <- NULL
    print(unlist(p))
    cat("\n")

    for (i in 1:length(object$secondary_models)) {
      cat(paste("Secondary model for parameter '", object$secondary_models[[i]]$par, "': ",
                object$secondary_models[[i]]$model, "\n", sep = ""))

      p <- object$secondary_models[[i]]
      p$model <- NULL
      p$par <- NULL

      cat(paste("Value under reference conditions: ", p$ref, "\n",
                sep = ""))
      p$ref <- NULL

      for (each_factor in names(p)) {

        cat(paste("Parameters with respect to ", each_factor, ":\n",
                  sep = ""))
        print(p[[each_factor]])
        cat("\n")

      }

      cat("\n")
    }

    cat(paste0("Treatment time: ", max(object$simulation$time, na.rm=TRUE), "\n"))
    cat(paste0("Minimum population size: ", min(object$simulation$logN, na.rm=TRUE), "\n"))
    cat(paste0("Total number of reductions: ",
               max(object$simulation$logN, na.rm=TRUE) - min(object$simulation$logN, na.rm=TRUE),
               "\n"
               ))

  }

}






















