
#' Dynamic treatment conditions
#' 
#' An example dataset to illustrate [predict_inactivation()] and [fit_inactivation()] 
#' under dynamic conditions. It describes a dummy dynamic treatment dependent on 
#' temperature and pH
#' 
#' @format A tibble with 3 rows and 3 columns:
#' \describe{
#'     \item{time}{elapsed time in min}
#'     \item{temp}{treatment temperature (C)}
#'     \item{pH}{treatment pH}
#' }
#' 
"dynamic_conditions"

#' Results of a dynamic treatment
#' 
#' An example dataset to illustrate [fit_inactivation()] 
#' under dynamic conditions. It describes simulated data observed under the 
#' dynamic treatment [dynamic_conditions].
#' 
#' @format A tibble with 20 rows and 3 columns:
#' \describe{
#'     \item{time}{elapsed time in min}
#'     \item{logN}{observed microbial concentration (log CFU/mL)}
#' }
#' 
"dynamic_experiment"

#' Results of an isothermal treatment
#' 
#' An example dataset to illustrate [fit_inactivation()] 
#' under isothermal conditions using simulated data
#' 
#' @format A tibble with 11 rows and 3 columns:
#' \describe{
#'     \item{time}{elapsed time in min}
#'     \item{logN}{observed microbial concentration (log CFU/mL)}
#' }
#' 
"inactivation_experiment"

#' Data for one-step fitting
#' 
#' An example dataset to illustrate [fit_inactivation()] 
#' using the one-step fitting approach. It simulates a batch of experiments at constant
#' environmental conditions changing two environmental factors.
#' 
#' @format A tibble with 25 rows and 4 columns:
#' \describe{
#'     \item{time}{treatment time in min}
#'     \item{logN}{observed microbial concentration (log CFU/mL)}
#'     \item{temp}{treatment temperature (C)}
#'     \item{pH}{treatment pH}
#' }
#' 
"inactivation_onestep"

#' Data for fitting a secondary inactivation model
#' 
#' An example dataset to illustrate [fit_inactivation_secondary()]. 
#' It simulates the value of `delta` estimated from a batch of experiments at 
#' constant conditions changing two environmental factors.
#' 
#' @format A tibble with 25 rows and 3 columns:
#' \describe{
#'     \item{temp}{treatment temperature (C)}
#'     \item{pH}{treatment pH}
#'     \item{delta}{the estimate of delta (min)}
#' }
#' 
"inactivation_secondary"

#' Data for two-step fitting
#' 
#' An example dataset to illustrate [fit_inactivation()] using the two-steps approach. 
#' It simulates the change in the microbial concentration observed in 9 independent experiments
#' at constant conditions changing two environmental factors.
#' 
#' @format A data.frame with 135 rows and 4 columns:
#' \describe{
#'     \item{temp}{treatment temperature (C)}
#'     \item{pH}{treatment pH}
#'     \item{time}{treatment time in min}
#'     \item{logN}{observed microbial concentration (log CFU/mL)}
#' }
#' 
"inactivation_twosteps"

#' Data for global fitting
#' 
#' An example dataset to illustrate [fit_inactivation()] using the global approach. 
#' It simulates the change in the microbial concentration observed in 2 independent experiments
#' at dynamic conditions. It should be used in combination with [multiple_environments]
#' 
#' @format A list of length two, with each element containing a tibble of 20 rows and 2 columns:
#' \describe{
#'     \item{time}{treatment time in min}
#'     \item{logN}{observed microbial concentration (log CFU/mL)}
#' }
#' 
"multiple_counts"

#' Environmental conditions for global fitting
#' 
#' An example dataset to illustrate [fit_inactivation()] using the global approach. 
#' It simulates the change in the environmental conditions in  2 independent experiments
#' at dynamic conditions. It should be used in combination with [multiple_counts]
#' 
#' @format A list of length two, with each element containing a tibble of 3 rows and 3 columns:
#' \describe{
#'     \item{time}{treatment time in min}
#'     \item{temp}{treatment temperature (C)}
#'     \item{pH}{treatment pH}
#' }
#' 
"multiple_environments"

