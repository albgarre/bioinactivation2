
#'
#'
sec_Bigelow <- function(x, xref, z) {
  - (x - xref)/z
}

#'
#'
sec_genBigelow <- function(x, xref, z, n) {
  - ((x - xref)/z)^n
}

#'
#'
sec_Lineal <- function(x, xref, b) {
  b*(x - xref)
}

#'
#'
sec_Arhenius <- function(x, xref, Ea) {
  R <- 8.31
  exp( Ea/R*(1/xref - 1/x) )
}

#'
#'
sec_logExponential <- function(x, k, xc) {

  log( 1 + exp( k*(x - xc) ) )
}
