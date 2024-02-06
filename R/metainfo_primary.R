
#' Metainformation of primary inactivation models
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Provides different types of meta-data about the primary inactivation models included
#' in bioinactivation2. This information is the basis of the automatic checks, and can also
#' help in the definition of models for [predict_inactivation()] and [fit_inactivation()].
#'
#' @param model_name The name of the model or `NULL` (default).
#'
#' @return
#' If model_name is `NULL`, returns a character string with the available models.
#' If is a valid identifier, it returns a list with meta-information about the model.
#' If model_name name is not a valid identifier, raises an error.
#'
#' @export
#'
primary_model_data <- function(model_name=NULL) {
  
  model_data <- list(
    Bigelow = list(identifier = "Bigelow",
                   name = "Bigelow (log-lineal) model",
                   pars = c("N0", "D"),
                   model = iso_Bigelow,
                   ref = paste(
                     "Bigelow, W. D. (1921). The Logarithmic Nature of Thermal Death Time Curves.",
                     "The Journal of Infectious Diseases, 29(5), 528-536."
                   )
    ),
    Peleg = list(
      identifier = "Peleg",
      name = "Peleg (Weibullian) model",
      pars = c("N0", "b", "n"),
      model = iso_Peleg,
      ref = paste(
        "Peleg, M., & Cole, M. B. (1998). Reinterpretation of Microbial Survival Curves.",
        "Critical Reviews in Food Science and Nutrition, 38(5), 353-380. https://doi.org/10.1080/10408699891274246"
      )
    ),
    Mafart = list(
      identifier = "Mafart",
      name = "Mafart (Weibullian) model",
      pars = c("N0", "delta", "p"),
      model = iso_Mafart,
      ref = paste(
        "Mafart, P., Couvert, O., Gaillard, S., & Leguerinel, I. (2002).",
        "On calculating sterility in thermal preservation methods: Application of the Weibull",
        "frequency distribution model. International Journal of Food Microbiology,",
        "72(1-2), 107-113. https://doi.org/10.1016/S0168-1605(01)00624-9"
      )
    ),
    Geeraerd = list(
      identifier = "Geeraerd",
      name = "Geeraerd model parameterized using the D-value",
      pars = c("N0", "D", "Nres", "SL"),
      model = iso_Geeraerd,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_k = list(
      identifier = "Geeraerd_k",
      name = "Geeraerd model parameterized using the inactivation rate (k)",
      pars = c("N0", "k", "Nres", "SL"),
      model = iso_Geeraerd_k,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_k_noTail = list(
      identifier = "Geeraerd_k_noTail",
      name = "Geeraerd model without tail parameterized using the inactivation rate (k)",
      pars = c("N0", "k", "SL"),
      model = iso_Geeraerd_noTail_k,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_noTail = list(
      identifier = "Geeraerd_noTail",
      name = "Geeraerd model without tail parameterized using the D-value",
      pars = c("N0", "D", "Nres"),
      model = iso_Geeraerd_noTail,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_k_noShoulder = list(
      identifier = "Geeraerd_k_noShoulder",
      name = "Geeraerd model without shoulder parameterized using the inactivation rate (k)",
      pars = c("N0", "k", "Nres"),
      model = iso_Geeraerd_noShoulder_k,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_noShoulder = list(
      identifier = "Geeraerd_noShoulder",
      name = "Geeraerd model without shoulder parameterized using the D-value",
      pars = c("N0", "D", "Nres"),
      model = iso_Geeraerd_noShoulder,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Metselaar = list(
      identifier = "Metselaar",
      name = "Metselaar (Weibullian) model",
      pars = c("N0", "D", "p", "Delta"),
      model = iso_Metselaar,
      ref = paste(
        "Metselaar, K. I., den Besten, H. M. W., Abee, T., Moezelaar, R., & Zwietering, M. H. (2013).",
        "Isolation and quantification of highly acid resistant variants of Listeria monocytogenes.",
        "International Journal of Food Microbiology, 166(3), 508-514. https://doi.org/10.1016/j.ijfoodmicro.2013.08.011"
      )
    ),
    Weibull_2phase = list(
      identifier = "Weibull_2phase",
      name = "2-population Weibullian inactivation model",
      pars = c("N0", "delta1", "delta2", "p1", "p2", "alpha"),
      model = iso_Weibull_2phase,
      ref = paste("")
    ),
    Trilinear = list(
      identifier = "Trilinear",
      name = "Trilineal inactivation model",
      pars = c("N0", "SL", "D", "Nres"),
      model = iso_Trilinear,
      ref = paste("")
    )
    
  )
  
  
  if (is.null(model_name)) {
    return(names(model_data))
  }
  
  my_model <- model_data[[model_name]]
  
  if (is.null(my_model)) {
    stop(paste("Unknown model name:", model_name))
  } else {
    my_model
  }
  
}


#' Basic check of parameters for primary models
#'
#' Checks that: the model name is correct, the right number of model
#' parameters have been defined and that the parameters have the right names
#'
#' @param model_name Model identifier
#' @param pars A named list of model parameters
#'
#' @return If there is no error, the model function.
#'
check_primary_pars <- function(model_name, pars) {
  
  ## (Indirectly) check that model name is correct
  
  my_data <- primary_model_data(model_name)
  model_pars <- my_data$pars
  
  ## Check the number of parameters
  
  if (length(model_pars) != length(pars)) {
    
    warning(paste0("The length of the parameters (", length(pars),
                   ") does not match the one of the model (", length(model_pars),
                   ").")
    )
  }
  
  ## Undo log-transformations 
  
  par_names <- str_replace(names(pars), "log", "")
  
  ## Check parameter names
  
  for(each_name in par_names) {
    
    if (!(each_name %in% model_pars)) {
      warning(paste("Not recognized parameter name:", each_name))
    }
    
  }
  
  ## Return
  
  my_data$model
  
}


