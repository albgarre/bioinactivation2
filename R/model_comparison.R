
#' Model comparison and selection for inactivation models
#' 
#' @description 
#' `r lifecycle::badge("experimental")`
#' 
#' This function is a constructor for [InactivationComparison] or [GlobalInactivationComparison],
#' a class that provides several functions for model comparison and model selection
#' for inactivation models fitted using [fit_inactivation()]. Please see the help pages for 
#' [InactivationComparison] or [GlobalInactivationComparison] for further details.
#' 
#' Although it is not necessary, we recommend passing the models as a named list,
#' as these names will later be kept in plots and tables.
#' 
#' @param models a (we recommend named) list of models fitted using [fit_inactivation()]. 
#' Every model should be of the same class. Otherwise, some functions may give unexpected results.
#' 
#' @importFrom FME modCost
#' 
#' @export
#' 
#' 
compare_inactivation_fits <- function(models) {
  
  # browser()
  
  ## Check for model types
  
  model_type <- unique(map_chr(models, ~ class(.)[1]))
  
  if (length(model_type) > 1) {
    warning("Every model should be of the same class. You are entering untested territory...")
  }
  
  ## Stop it if this is a two-step fitting
  
  if (models[[1]]$approach == "two-steps") {
    stop("Model comparison for two-steps fitting not implemented")
  }

  ## Calculate residuals

  if (models[[1]]$approach == "global") {
    
    residuals <- models %>%
      map(
        ~ imap_dfr(.$data, ~ mutate(.x, exp = .y))
      ) %>%
      map2(models, ~ mutate(.x, res = residuals(.y)))
    
  } else if (models[[1]]$approach == "one-step") {
    
    residuals <- models %>%
      map(
        ~ mutate(.$data, pred = fitted(.), res = residuals(.))
      ) %>%
      imap(~ mutate(.x, model = .y))

  } else {

    d <- as.data.frame(models[[1]]$data)

    t <- seq(0, max(d$time, na.rm = TRUE), length = 1000)

    residuals <- models %>%
      map(
        ~ data.frame(time = t,
                     logN = predict(., times = t)
        )
      ) %>%
      map(~ modCost(model = ., obs = d)
      )

  }

  ## Save the type of environment

  environment <- models[[1]]$environment

  ## Save the type of algorithm

  algorithm <- models[[1]]$algorithm

  ## Return

  out <- list(models = models,
              residuals = residuals,
              environment = environment,
              algorithm = algorithm,
              approach = models[[1]]$approach)

  if (models[[1]]$approach == "global") {

    class(out) <- c("GlobalComparison", class(out))

  } else if (models[[1]]$approach == "one-step") {
    
    class(out) <- c("OneStepComparison", class(out))
    
  } else {

    class(out) <- c("InactivationComparison", class(out))

  }

  out
  
}