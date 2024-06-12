
#' Metainformation of secondary inactivation models
#' 
#' @description 
#' `r lifecycle::badge("stable")`
#' 
#' Provides different types of meta-data about the secondary inactivation models included
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
secondary_model_data <- function(model_name=NULL) {
  
  model_data <- list(
    Bigelow = list(identifier = "Bigelow",
                   name = "Bigelow secondary model",
                   pars = c("xref", "z"),
                   model = sec_Bigelow,
                   ref = paste(
                     "Bigelow, W. D., & Esty, J. R. (1920).",
                     "The thermal death point in relation to time of typical thermophilic organisms.",
                     "The Journal of Infectious Diseases, 602-617."
                   )
    ),
    genBigelow = list(identifier = "genBigelow",
                      name = "Generalized-Bigelow secondary model",
                      pars = c("xref", "z", "n"),
                      model = sec_genBigelow,
                      ref = paste("")
    ),
    Lineal = list(identifier = "Lineal",
                  name = "Lineal secondary model",
                  pars = c("xref", "b"),
                  model = sec_Lineal,
                  ref = paste("")
    ),
    Arrhenius = list(identifier = "Arrhenius",
                     name = "Arrhenius-type secondary model",
                     pars = c("xref", "Ea"),
                     model = sec_Arrhenius,
                     ref = paste("")
    ),
    logExponential = list(identifier = "logExponential",
                          name = "log-exponential secondary model",
                          pars = c("k", "xc"),
                          model = sec_logExponential,
                          ref = paste(
                            "Peleg, M., & Cole, M. B. (1998). Reinterpretation of Microbial Survival Curves.",
                            "Critical Reviews in Food Science and Nutrition, 38(5), 353-380. https://doi.org/10.1080/10408699891274246"
                          )
    ),
    kAcclimation = list(identifier = "kAcclimation",
                          name = "secondary model for the acclimation model",
                          pars = c("E", "Xsi"),
                          model = sec_Acclimation,
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

#' Basic check of parameters for fitting secondary models
#'
#' Checks that: the model name is correct, the right number of model
#' parameters have been defined and that the parameters have the right names
#'
#' @param model_name model identifier
#' @param pars a named list of model parameters
#' @param vars a character vector of environmental factors included in the model
#'
#' @return If there is no error, the model function.
#'
check_secondary_pars <- function(model_name, pars, vars) {
  
  ## (Indirectly) check that model name is correct
  
  my_data <- secondary_model_data(model_name)
  model_pars <- my_data$pars
  
  ## Check the reference conditions
  
  if (!("ref" %in% names(pars) || "logref" %in% names(pars))) {
    warning("Value at reference conditions not identified")
  }
  
  ## Do checks for each factor
  
  lapply(vars, function(each_factor) {
    
    ## Get the parameter for this factor
    
    this_p <- pars[grepl(paste0(each_factor, "_"), names(pars))]
    
    ## Check the number of parameters
    
    if (length(model_pars) != length(this_p)) {
      
      warning(paste0("The length of the parameters (", length(this_p),
                     ") does not match the one of the model (", length(model_pars),
                     ").")
      )
    }
    
    ## Remove the factor name from the parameter
    
    par_names <- str_replace(names(this_p), paste0(each_factor, "_"), "")
    
    ## Undo log-transformations 
    
    par_names <- str_replace(par_names, "log", "")
    
    ## Check parameter names
    
    for(each_name in par_names) {
      
      if (!(each_name %in% model_pars)) {
        warning(paste("Not recognized parameter name:", each_name))
      }
      
    }
    
    
  })
  
  ## Return
  
  my_data$model
  
}





