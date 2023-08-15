
#'
#'
iso_Bigelow <- function(times, N0, D) {
  logN0 <- log10(N0)
  logN0 - times/D
}

#'
#'
iso_Peleg <- function(times, N0, b, n) {
  logN0 <- log10(N0)
  logN0 - b*times^n
}

#'
#'
iso_Mafart <- function(times, N0, delta, p) {
  logN0 <- log10(N0)
  logN0 - (times/delta)^p
}

#'
#'
iso_Geeraerd = function(times, N0, D, Nres, SL) {

  k <- log(10)/D

  iso_Geeraerd_k(times, N0, SL, k, Nres)

}

#'
#'
iso_Metselaar = function(times, Delta, N0, D, p) {

  logN0 <- log10(N0)

  p$logN0 - Delta*(times/Delta/D)^p

}

#'
#'
iso_Weibull_2phase = function(times, N0, delta1, delta2, p1, p2, alpha) {

  part1 <- 10^(- (times/delta1)^p1 + alpha )
  part2 <- 10^(- (times/delta2)^p2)
  N <- N0/(1 + 10^alpha)*(part1 + part2)

  log10(N)

}

#'
#'
iso_Trilinear = function(times, N0, SL, D, Nres) {

  logN0 <- log10(N0)
  logNres <- log10(Nres)

  tibble(logN = logN0 - (times - SL)/D) %>%
    mutate(logN = ifelse(times < SL, logN0, logN),
           logN = ifelse(logN < logNres, logNres, logN)) %>%
    pull(logN)
}

#'
#'
iso_Geeraerd_k = function(times, N0, SL, k, Nres) {

  logN0 <- log10(N0)
  logNres <- log10(Nres)

  logNres + log10(( (10^(logN0-logNres)-1)*exp(k*SL) )/(exp(k*times) + exp(k*SL) - 1) + 1)

}

iso_Geeraerd_noTail_k = function(times, N0, SL, k) {

  logN0 <- log10(N0)

  N <- 10^logN0 * exp(-k*times) * exp(k*SL) / ( 1 + ( exp(k*SL) - 1 )*exp(-k*times) )

  log10(N)

}

#'
#'
iso_Geeraerd_noTail = function(times, N0, SL, D) {

  k <- log(10)/D

  iso_Geeraerd_noTail_k(times, N0, SL, k)

}

#'
#'
iso_Geeraerd_noShoulder_k = function(times, N0, k, Nres) {

  logN0 <- log10(N0)
  logNres <- log10(Nres)

  logNres + log10(( (10^(logN0 - logNres) - 1))/(exp(k*times)) + 1)

}

iso_Geeraerd_noShoulder = function(times, N0, D, Nres) {

  k <- log(10)/D

  iso_Geeraerd_noShoulder_k(times, N0, k, Nres)

}
