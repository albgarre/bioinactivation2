
#' Metainformation of primary inactivation models for dynamic conditions
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Provides different types of meta-data about the primary inactivation models included
#' in bioinactivation2 for dynamic calculations. This information is the basis of the automatic checks, and can also
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
dynamic_model_data <- function(model_name=NULL) {
  
  model_data <- list(
    Bigelow = list(identifier = "Bigelow",
                   name = "Bigelow (1st order kinetics) model",
                   pars_primary = c("N0"),
                   pars_secondary = c("D"),
                   model = dyna_Bigelow,
                   ref = paste(
                     "Bigelow, W. D. (1921). The Logarithmic Nature of Thermal Death Time Curves.",
                     "The Journal of Infectious Diseases, 29(5), 528-536."
                   )
    ),
    Peleg = list(
      identifier = "Peleg",
      name = "Peleg model",
      pars_primary = c("N0"),
      pars_secondary = c("b", "n"),
      model = dyna_Peleg,
      ref = paste(
        "Peleg, M., & Cole, M. B. (1998). Reinterpretation of Microbial Survival Curves.",
        "Critical Reviews in Food Science and Nutrition, 38(5), 353-380. https://doi.org/10.1080/10408699891274246"
      )
    ),
    Mafart = list(
      identifier = "Mafart",
      name = "Mafart model",
      pars_primary = c("N0"),
      pars_secondary = c("delta", "p"),
      model = dyna_Mafart,
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
      pars_primary = c("N0", "C0"),
      pars_secondary = c("D", "Nres"),
      model = dyna_Geeraerd,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_k = list(
      identifier = "Geeraerd_k",
      name = "Geeraerd model parameterized using the inactivation rate (k)",
      pars_primary = c("N0", "C0"),
      pars_secondary = c("k", "Nres"),
      model = dyna_Geeraerd_k,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_k_noTail = list(
      identifier = "Geeraerd_k_noTail",
      name = "Geeraerd model without tail parameterized using the inactivation rate (k)",
      pars_primary = c("N0", "C0"),
      pars_secondary = c("k"),
      model = dyna_Geeraerd_noTail_k,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_noTail = list(
      identifier = "Geeraerd_noTail",
      name = "Geeraerd model without tail parameterized using the D-value",
      pars_primary = c("N0", "C0"),
      pars_secondary = c("D"),
      model = dyna_Geeraerd_noTail,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_k_noShoulder = list(
      identifier = "Geeraerd_k_noShoulder",
      name = "Geeraerd model without shoulder parameterized using the inactivation rate (k)",
      pars_primary = c("N0"),
      pars_secondary = c("k", "Nres"),
      model = dyna_Geeraerd_noSL_k,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Geeraerd_noShoulder = list(
      identifier = "Geeraerd_noShoulder",
      name = "Geeraerd model without shoulder parameterized using the D-value",
      pars_primary = c("N0"),
      pars_secondary = c("D", "Nres"),
      model = dyna_Geeraerd_noSL,
      ref = paste(
        "Geeraerd, A. H., Herremans, C. H., & Van Impe, J. F. (2000).",
        "Structural model requirements to describe microbial inactivation during a mild heat treatment.",
        "International Journal of Food Microbiology, 59(3), 185-209. https://doi.org/10.1016/S0168-1605(00)00362-7"
      )
    ),
    Acclimation = list(
      identifier = "Acclimation",
      name = "Modification of the Bigelow model considering stress acclimation",
      pars_primary = c("N0", "p0"),
      pars_secondary = c("D", "c", "k"),
      model = dyna_acclimation,
      ref = paste(
        "Garre, A., Huertas, J. P., González-Tejedor, G. A., Fernández, P. S., Egea, J. A., Palop, A., & Esnoz, A. (2018).",
        "Mathematical quantification of the induced stress resistance of microbial populations during",
        "non-isothermal stresses. International Journal of Food Microbiology,",
        "266, 133–141. https://doi.org/10.1016/j.ijfoodmicro.2017.11.023"
        
      )
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


#' Basic check of parameters for predictions under dynamic conditions
#'
#' Checks that: AA
#'
#' @param primary_model A named list as in [predict_inactivation()]
#' @param secondary_models A named list as in [predict_inactivation()]
#' 
#' @importFrom purrr %>% map_chr
#'
#' @return If there is no error, the model function.
#'
check_dynamic_pars <- function(primary_model, secondary_models) {
  
  # browser()
  
  ## (Indirectly) check that model name is correct
  
  my_data <- dynamic_model_data(primary_model$model)
  
  ## Check parameters of the primary model
  
  primary_model$model <- NULL
  par_names <- names(primary_model)
  par_names <- str_replace(par_names, "log", "")  # Undo log-transformations
  
  model_pars <- my_data$pars_primary
  
  #- Check the number of parameters
  
  if (length(model_pars) != length(par_names)) {
    
    warning(paste0("The length of the parameters (", length(par_names),
                   ") does not match the one of the model (", length(model_pars),
                   ").")
    )
  }
  
  #- Check parameter names
  
  for(each_name in par_names) {
    
    if (!(each_name %in% model_pars)) {
      warning(paste("Not recognized parameter name:", each_name))
    }
    
  }
  
  ## Check every parameter has a secondary model
  
  model_pars <- my_data$pars_secondary
  
  par_names <- secondary_models %>% map_chr(~.$par)
  par_names <- str_replace(par_names, "log", "")  # Undo log-transformations
  
  #- Check the number of parameters
  
  if (length(model_pars) != length(par_names)) {
    
    warning(paste0("The length of the parameters (", length(par_names),
                   ") does not match the one of the model (", length(model_pars),
                   ").")
    )
  }
  
  #- Check parameter names
  
  for(each_name in par_names) {
    
    if (!(each_name %in% model_pars)) {
      warning(paste("Not recognized parameter name:", each_name))
    }
    
  }
  
  ## Check each secondary model
  
  lapply(secondary_models, function(this_model) {
    
    # browser()
    
    ## (Indirectly) check that model name is correct
    
    my_data <- secondary_model_data(this_model$model)
    model_pars <- my_data$pars
    
    ## Check the reference conditions
    
    if ( !("ref" %in% names(this_model)) ) {
      warning(paste("Value at reference conditions not identified for parameter:", 
                    this_model$par))
    }
    
    ## Do checks for each factor
    
    env_models <- this_model
    env_models$par <- NULL
    env_models$model <- NULL
    env_models$ref <- NULL
    
    lapply(env_models, function(this_factor) {
      
      ## Get the parameter for this factor
      
      this_p <- names(this_factor)
      
      ## Check the number of parameters
      
      if (length(model_pars) != length(this_p)) {
        
        warning(paste0("The length of the parameters (", length(this_p),
                       ") does not match the one of the model (", length(model_pars),
                       ").")
        )
      }
      
      ## Undo log-transformations 
      
      this_p <- str_replace(this_p, "log", "")
      
      ## Check parameter names
      
      for(each_name in this_p) {
        
        if (!(each_name %in% model_pars)) {
          warning(paste("Not recognized parameter name:", each_name))
        }
        
      }
      
      
    })
      
  })
  
  ## Return
  
  my_data$model
  
}




