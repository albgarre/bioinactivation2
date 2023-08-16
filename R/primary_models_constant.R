
#' Bigelow (log-lineal) model for constant conditions
#' 
#' @param times a numerical vector of time points
#' @param N0 initial microbial concentration (1D numeric)
#' @param D D-value (1D numeric)
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
#'
iso_Bigelow <- function(times, N0, D) {
  logN0 <- log10(N0)
  logN0 - times/D
}

#' Peleg (Weibullian) model for constant conditions
#' 
#' @param times a numerical vector of time points
#' @param N0 initial microbial concentration (1D numeric)
#' @param b b-value of the Peleg model (1D numeric)
#' @param n shape factor of the Peleg model (1D numeric)
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
iso_Peleg <- function(times, N0, b, n) {
  logN0 <- log10(N0)
  logN0 - b*times^n
}

#' Mafart (Weibullian) model for constant conditions
#' 
#' @param times a numerical vector of time points
#' @param N0 initial microbial concentration (1D numeric)
#' @param delta delta-value (1D numeric)
#' @param p shape factor of the Mafart model (1D numeric)
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
iso_Mafart <- function(times, N0, delta, p) {
  logN0 <- log10(N0)
  logN0 - (times/delta)^p
}

#' Geeraerd model for constant conditions (based on the D-value)
#' 
#' @param times a numerical vector of time points
#' @param N0 initial microbial concentration (1D numeric)
#' @param D D-value (1D numeric)
#' @param Nres microbial concentration at the tail (1D numeric)
#' @param SL shoulder length (1D numeric)
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
iso_Geeraerd = function(times, N0, D, Nres, SL) {

  k <- log(10)/D

  iso_Geeraerd_k(times, N0, SL, k, Nres)

}

#' Metselaar (Weibull) model for constant conditions
#' 
#' @param times a numerical vector of time points
#' @param N0 initial microbial concentration (1D numeric)
#' @param Delta Reference log-reductions for the model (1D numeric)
#' @param D Equivalent D-value for Delta log-reductions (1D numeric)
#' @param p shape factor of the Weibull model (1D numeric)
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
iso_Metselaar = function(times, Delta, N0, D, p) {

  logN0 <- log10(N0)

  p$logN0 - Delta*(times/Delta/D)^p

}

#' Weibull model with two populations
#' 
#' @param times a numerical vector of time points
#' @param N0 initial microbial concentration (1D numeric)
#' @param delta1 delta-value of the 1st population (1D numeric)
#' @param delta2 delta-value of the 2nd population (1D numeric)
#' @param p1 shape factor of the Weibull model for the 1st population (1D numeric)
#' @param p2 shape factor of the Weibull model for the 2nd population (1D numeric)
#' @param alpha parameter describing the population ratio (1D numeric)
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
iso_Weibull_2phase = function(times, N0, delta1, delta2, p1, p2, alpha) {

  part1 <- 10^(- (times/delta1)^p1 + alpha )
  part2 <- 10^(- (times/delta2)^p2)
  N <- N0/(1 + 10^alpha)*(part1 + part2)

  log10(N)

}

#' Trilineal inactivatino model
#' 
#' @param times a numerical vector of time points
#' @param N0 initial microbial concentration (1D numeric)
#' @param D D-value (1D numeric)
#' @param Nres microbial concentration at the tail (1D numeric)
#' @param SL shoulder length (1D numeric)
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
iso_Trilinear = function(times, N0, SL, D, Nres) {

  logN0 <- log10(N0)
  logNres <- log10(Nres)

  tibble(logN = logN0 - (times - SL)/D) %>%
    mutate(logN = ifelse(times < SL, logN0, logN),
           logN = ifelse(logN < logNres, logNres, logN)) %>%
    pull(logN)
}

#' Geeraerd model for constant conditions (based on k)
#' 
#' @param times a numerical vector of time points
#' @param N0 initial microbial concentration (1D numeric)
#' @param k inactivation rate (1D numeric)
#' @param Nres microbial concentration at the tail (1D numeric)
#' @param SL shoulder length (1D numeric)
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
iso_Geeraerd_k = function(times, N0, SL, k, Nres) {

  logN0 <- log10(N0)
  logNres <- log10(Nres)

  logNres + log10(( (10^(logN0-logNres)-1)*exp(k*SL) )/(exp(k*times) + exp(k*SL) - 1) + 1)

}

#' Geeraerd model for constant conditions without tail (based on k)
#' 
#' @inheritParams iso_Geeraerd_k
#' 
#' @returns A numerical vector of log-microbial concentrations
#' 
iso_Geeraerd_noTail_k = function(times, N0, SL, k) {

  logN0 <- log10(N0)

  N <- 10^logN0 * exp(-k*times) * exp(k*SL) / ( 1 + ( exp(k*SL) - 1 )*exp(-k*times) )

  log10(N)

}

#' Geeraerd model for constant conditions without tail (based on K)
#' 
#' @inheritParams iso_Geeraerd
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
iso_Geeraerd_noTail = function(times, N0, SL, D) {

  k <- log(10)/D

  iso_Geeraerd_noTail_k(times, N0, SL, k)

}

#' Geeraerd model for constant conditions without shoulder (based on k)
#' 
#' @inheritParams iso_Geeraerd_k
#' 
#' @returns A numerical vector of log-microbial concentrations
#'
iso_Geeraerd_noShoulder_k = function(times, N0, k, Nres) {

  logN0 <- log10(N0)
  logNres <- log10(Nres)

  logNres + log10(( (10^(logN0 - logNres) - 1))/(exp(k*times)) + 1)

}

#' Geeraerd model for constant conditions without tail (based on D)
#' 
#' @inheritParams iso_Geeraerd
#' 
#' @returns A numerical vector of log-microbial concentrations
#' 
iso_Geeraerd_noShoulder = function(times, N0, D, Nres) {

  k <- log(10)/D

  iso_Geeraerd_noShoulder_k(times, N0, k, Nres)

}




