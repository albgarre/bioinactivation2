
#' Bigelow-type secondary model
#' 
#' @param x value of the environmental condition
#' @param xref reference value of x
#' @param z z-value
#'
sec_Bigelow <- function(x, xref, z) {
  - (x - xref)/z
}


#' Genearlized-Bigelow secondary model
#' 
#' @param x value of the environmental condition
#' @param xref reference value of x
#' @param z z-value
#' @param n order of the model
#'
#'
sec_genBigelow <- function(x, xref, z, n) {
  - ((x - xref)/z)^n
}


#' Lineal secondary model
#' 
#' @param x value of the environmental condition
#' @param xref reference value of x
#' @param b slope of the relation
#'
#'
sec_Lineal <- function(x, xref, b) {
  b*(x - xref)
}


#' Arrhenius-type secondary model
#' 
#' @param x value of the environmental condition
#' @param xref reference value of x
#' @param Ea Activation energy
#'
#'
sec_Arrhenius <- function(x, xref, Ea) {
  R <- 8.31
  exp( Ea/R*(1/xref - 1/x) )
}


#' Log-exponential secondary model
#' 
#' @param x value of the environmental condition
#' @param xc critical value of x
#' @param k k-parameter of the secondary model
#'
#'
sec_logExponential <- function(x, k, xc) {

  log( 1 + exp( k*(x - xc) ) )
}






