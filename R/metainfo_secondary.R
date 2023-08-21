
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
                     "The Journal of Infectious Diseases, 602–617."
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
                     pars = c("xref", "Ea"),
                     model = sec_logExponential,
                     ref = paste(
                       "Peleg, M., & Cole, M. B. (1998). Reinterpretation of Microbial Survival Curves.",
                       "Critical Reviews in Food Science and Nutrition, 38(5), 353–380. https://doi.org/10.1080/10408699891274246"
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