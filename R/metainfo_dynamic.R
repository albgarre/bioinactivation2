
#' Metainformation of primary inactivation models for dynamic conditions
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
dynamic_model_data <- function(model_name=NULL) {
  
  # model_data <- list(
  #   Bigelow = list(identifier = "Bigelow",
  #                  name = "Bigelow (log-lineal) model",
  #                  pars = c("N0", "D"),
  #                  model = iso_Bigelow,
  #                  ref = paste(
  #                    "Bigelow, W. D. (1921). The Logarithmic Nature of Thermal Death Time Curves.",
  #                    "The Journal of Infectious Diseases, 29(5), 528–536."
  #                  )
  #   )
  #   
  # )
  # 
  # 
  # if (is.null(model_name)) {
  #   return(names(model_data))
  # }
  # 
  # my_model <- model_data[[model_name]]
  # 
  # if (is.null(my_model)) {
  #   stop(paste("Unknown model name:", model_name))
  # } else {
  #   my_model
  # }
  
}




