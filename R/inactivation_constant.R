

#' Predicting inactivation under constant environmental conditions
#' 
#' Predictions using only primary inactivation models. In principle, these models
#' would be valid for constant environmental conditions.
#' 
#' @param times a numeric vector of time points where to calculate the prediction
#' @param model_name a model identifier as per [primary_model_data()]
#' @param model_pars a numeric vector (or list) describing the values of the model
#' parameters as defined by [primary_model_data()]
#' @param check whether to due some model checks. `TRUE` by default
#' 
#' @returns a tibble with two columnes ('time' and 'logN') with the results of the
#' calculation.
#'
predict_constant_inactivation <- function(times,
                                          model_name,
                                          model_pars,
                                          check = TRUE
                                          ) {
  
  ## Check the model parameters

  model_pars <- as.list(model_pars)

  if (isTRUE(check)) {

    check_primary_pars(model_name, model_pars)

  }

  ## Calculate the prediction

  logN <- switch(model_name,
                 Bigelow = iso_Bigelow(times,
                                       model_pars$N0,
                                       model_pars$D),
                 Peleg = iso_Peleg(times,
                                   model_pars$N0,
                                   model_pars$b,
                                   model_pars$n),
                 Mafart = iso_Mafart(times,
                                     model_pars$N0,
                                     model_pars$delta,
                                     model_pars$p
                 ),
                 Geeraerd = iso_Geeraerd(times,
                                         model_pars$N0,
                                         model_pars$D,
                                         model_pars$Nres,
                                         model_pars$SL
                 ),
                 Metselaar = iso_Metselaar(times,
                                           model_pars$Delta,
                                           model_pars$N0,
                                           model_pars$D,
                                           model_pars$p),
                 Weibull_2phase = iso_Weibull_2phase(times,
                                                     model_pars$N0,
                                                     model_pars$delta1,
                                                     model_pars$delta2,
                                                     model_pars$p1,
                                                     model_pars$p2,
                                                     model_pars$alpha),
                 Trilinear = iso_Trilinear(times,
                                           model_pars$N0,
                                           model_pars$SL,
                                           model_pars$D,
                                           model_pars$Nres),
                 Geeraerd_k = iso_Geeraerd_k(times,
                                             model_pars$N0,
                                             model_pars$SL,
                                             model_pars$k,
                                             model_pars$Nres),
                 Geeraerd_k_noTail = iso_Geeraerd_noTail_k(times,
                                                           model_pars$N0,
                                                           model_pars$SL,
                                                           model_pars$k),
                 Geeraerd_noTail = iso_Geeraerd_noTail(times,
                                                       model_pars$N0,
                                                       model_pars$SL,
                                                       model_pars$D),
                 Geeraerd_k_noShoulder = iso_Geeraerd_noShoulder_k(times,
                                                                   model_pars$N0,
                                                                   model_pars$k,
                                                                   model_pars$Nres),
                 Geeraerd_noShoulder = iso_Geeraerd_noShoulder(times,
                                                               model_pars$N0,
                                                               model_pars$D,
                                                               model_pars$Nres),
                 stop(paste("Unknown model:", model_name))
  )

  ## Prepare the output

  my_sim <- tibble(time = times, logN = logN)

  my_sim

}














