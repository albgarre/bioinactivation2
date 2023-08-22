
#' Bigelow model for predictions under dynamic conditions
#' 
#' The functions is defined to be called from [ode()] from the deSolve package.
#' 
#' @param time a 1D numeric with the treatment time
#' @param state a named numeric vector with one element defining the microbial concentration
#' @param primary_pars ignored. Just a dumpster for the param argument from ode
#' @param env_interpolator  an interpolator for the environmental conditions as returned by [approx_env()]
#' @param secondary_models a nested list defining the secondary models as per [apply_secondary_models()]
#' 
#' @returns a list with dN/dt as per [ode()]
#'
dyna_Bigelow <- function(time, state, primary_pars, env_interpolator, secondary_models) {

  ## Get the environmental conditions

  env_conditions <-  map(env_interpolator, ~.(time))

  ## Apply the secondary models

  p <- apply_secondary_models(env_conditions, secondary_models)

  state <- as.list(state)

  dN <- -log(10)*state$N/p$D

  list(c(N = dN))

}


#' Mafart model for predictions under dynamic conditions
#' 
#' The functions is defined to be called from [ode()] from the deSolve package.
#' 
#' The Mafart model under dynamic conditions can be singular for t=0 for some combinations
#' of model parameters. This is fixed by adding a tolerance to every time point. By default,
#' the tolerance is 1e-12, so it should not affect the calculations. Nonetheless,
#' this can be changed using the `tol_time` argument.
#' 
#' @param time a 1D numeric with the treatment time
#' @param state a named numeric vector with one element defining the microbial concentration
#' @param primary_pars ignored. Just a dumpster for the param argument from ode
#' @param env_interpolator  an interpolator for the environmental conditions as returned by [approx_env()]
#' @param secondary_models a nested list defining the secondary models as per [apply_secondary_models()]
#' @param tol_time a tolerance for the time (to avoid singularities for t=0). See Details
#' 
#' @returns a list with dN/dt as per [ode()]
#'
dyna_Mafart <- function(time, state, primary_pars, env_interpolator, secondary_models,
                        tol_time = 1e-12) {

  ## Get the environmental conditions

  env_conditions <-  map(env_interpolator, ~.(time))

  ## Apply the secondary models

  p <- apply_secondary_models(env_conditions, secondary_models)

  ## Apply the tolerance to avoid numerical issues

  time <- time + tol_time

  ## Calculate

  state <- as.list(state)

  dN <- - state$N * p$p * (1/p$delta)^p$p * time^(p$p-1) * log(10)

  list(c(N = dN))

}

#' Peleg model for predictions under dynamic conditions
#' 
#' The functions is defined to be called from [ode()] from the deSolve package.
#' 
#' The Peleg model under dynamic conditions can be singular for t=0 for some combinations
#' of model parameters. This is fixed by adding a tolerance to every time point. By default,
#' the tolerance is 1e-12, so it should not affect the calculations. Nonetheless,
#' this can be changed using the `tol_time` argument.
#' 
#' @param time a 1D numeric with the treatment time
#' @param state a named numeric vector with one element defining the microbial concentration
#' @param primary_pars ignored. Just a dumpster for the param argument from ode
#' @param env_interpolator  an interpolator for the environmental conditions as returned by [approx_env()]
#' @param secondary_models a nested list defining the secondary models as per [apply_secondary_models()]
#' @param tol_time a tolerance for the time (to avoid singularities for t=0). See Details
#' 
#' @returns a list with dN/dt as per [ode()]
#'
#'
dyna_Peleg <- function(time, state, primary_pars, env_interpolator, secondary_models,
                       tol_logS = 1e-12) {

  ## Get the environmental conditions

  env_conditions <-  map(env_interpolator, ~.(time))

  ## Apply the secondary models

  p <- apply_secondary_models(env_conditions, secondary_models)

  # print(p)
  # print(env_conditions)

  state <- as.list(state)

  logS <- log10(state$N) - log10(primary_pars$N0)

  ## Apply the tolerance to avoid numerical issues

  logS <- logS - tol_logS

  ## Calculate the derivative

  dlogS <- -p$b * p$n * ( - logS/p$b ) ^( (p$n-1)/p$n)

  dN <- dlogS*state$N*log(10)

  list(c(N = dN))

}

#' Geeraerd model for predictions under dynamic conditions
#' 
#' The functions is defined to be called from [ode()] from the deSolve package.
#' 
#' @param time a 1D numeric with the treatment time
#' @param state a named numeric vector with two element defining the microbial concentration and 
#' the value of the ideal substance 'C'
#' @param primary_pars ignored. Just a dumpster for the param argument from ode
#' @param env_interpolator  an interpolator for the environmental conditions as returned by [approx_env()]
#' @param secondary_models a nested list defining the secondary models as per [apply_secondary_models()]
#' @param tol_time a tolerance for the time (to avoid singularities for t=0). See Details
#' 
#' @returns a list with two elements: dN/dt and dC/dt, as per [ode()]
#'
#'
dyna_Geeraerd <- function(time, state, primary_pars, env_interpolator, secondary_models) {
  
  
  # browser()
  
  ## Get the environmental conditions
  
  env_conditions <-  map(env_interpolator, ~.(time))
  
  ## Apply the secondary models
  
  p <- apply_secondary_models(env_conditions, secondary_models)
  
  state <- as.list(state)
  
  dC <- -log(10)*state$C/p$D
  
  alpha <- 1/(1 + state$C)
  beta <- (1 - p$Nres/state$N)
  
  dN <- -log(10)*state$N/p$D * alpha * beta
  
  res <- c(dN, dC)
  return(list(res))
  
}

#' Geeraerd model for predictions under dynamic conditions parameterized using k
#' 
#' The functions is defined to be called from [ode()] from the deSolve package.
#' 
#' @param time a 1D numeric with the treatment time
#' @param state a named numeric vector with two element defining the microbial concentration and 
#' the value of the ideal substance 'C'
#' @param primary_pars ignored. Just a dumpster for the param argument from ode
#' @param env_interpolator  an interpolator for the environmental conditions as returned by [approx_env()]
#' @param secondary_models a nested list defining the secondary models as per [apply_secondary_models()]
#' @param tol_time a tolerance for the time (to avoid singularities for t=0). See Details
#' 
#' @returns a list with two elements: dN/dt and dC/dt, as per [ode()]
#'
#'
dyna_Geeraerd_k <- function(time, state, primary_pars, env_interpolator, secondary_models) {
  
  ## Get the environmental conditions
  
  env_conditions <-  map(env_interpolator, ~.(time))
  
  ## Apply the secondary models
  
  p <- apply_secondary_models(env_conditions, secondary_models)
  
  state <- as.list(state)
  
  dC <- - p$k * state$C
  
  alpha <- 1/(1 + state$C)
  beta <- (1 - p$Nres/state$N)
  
  dN <- - p$k * state$N * alpha * beta
  
  res <- c(dN, dC)
  return(list(res))
  
}

#' Geeraerd model without tail for predictions under dynamic conditions
#' 
#' The functions is defined to be called from [ode()] from the deSolve package.
#' 
#' @param time a 1D numeric with the treatment time
#' @param state a named numeric vector with two element defining the microbial concentration and 
#' the value of the ideal substance 'C'
#' @param primary_pars ignored. Just a dumpster for the param argument from ode
#' @param env_interpolator  an interpolator for the environmental conditions as returned by [approx_env()]
#' @param secondary_models a nested list defining the secondary models as per [apply_secondary_models()]
#' @param tol_time a tolerance for the time (to avoid singularities for t=0). See Details
#' 
#' @returns a list with two elements: dN/dt and dC/dt, as per [ode()]
#'
#'
dyna_Geeraerd_noTail <- function(time, state, primary_pars, env_interpolator, secondary_models) {
  
  
  # browser()
  
  ## Get the environmental conditions
  
  env_conditions <-  map(env_interpolator, ~.(time))
  
  ## Apply the secondary models
  
  p <- apply_secondary_models(env_conditions, secondary_models)
  
  state <- as.list(state)
  
  dC <- -log(10)*state$C/p$D
  
  alpha <- 1/(1 + state$C)
  beta <- 1
  
  dN <- -log(10)*state$N/p$D * alpha * beta
  
  res <- c(dN, dC)
  return(list(res))
  
}

#' Geeraerd model without a tail for predictions under dynamic conditions parameterized on k
#' 
#' The functions is defined to be called from [ode()] from the deSolve package.
#' 
#' @param time a 1D numeric with the treatment time
#' @param state a named numeric vector with two element defining the microbial concentration and 
#' the value of the ideal substance 'C'
#' @param primary_pars ignored. Just a dumpster for the param argument from ode
#' @param env_interpolator  an interpolator for the environmental conditions as returned by [approx_env()]
#' @param secondary_models a nested list defining the secondary models as per [apply_secondary_models()]
#' @param tol_time a tolerance for the time (to avoid singularities for t=0). See Details
#' 
#' @returns a list with two elements: dN/dt and dC/dt, as per [ode()]
#'
#'
dyna_Geeraerd_noTail_k <- function(time, state, primary_pars, env_interpolator, secondary_models) {
  
  ## Get the environmental conditions
  
  env_conditions <-  map(env_interpolator, ~.(time))
  
  ## Apply the secondary models
  
  p <- apply_secondary_models(env_conditions, secondary_models)
  
  state <- as.list(state)
  
  dC <- - p$k * state$C
  
  alpha <- 1/(1 + state$C)
  beta <- 1
  
  dN <- - p$k * state$N * alpha * beta
  
  res <- c(dN, dC)
  return(list(res))
  
}

#' Geeraerd model without shoulder for predictions under dynamic conditions
#' 
#' The functions is defined to be called from [ode()] from the deSolve package.
#' 
#' @param time a 1D numeric with the treatment time
#' @param state a named numeric vector with two element defining the microbial concentration and 
#' the value of the ideal substance 'C'
#' @param primary_pars ignored. Just a dumpster for the param argument from ode
#' @param env_interpolator  an interpolator for the environmental conditions as returned by [approx_env()]
#' @param secondary_models a nested list defining the secondary models as per [apply_secondary_models()]
#' @param tol_time a tolerance for the time (to avoid singularities for t=0). See Details
#' 
#' @returns a list with two elements: dN/dt and dC/dt, as per [ode()]
#'
#'
dyna_Geeraerd_noSL <- function(time, state, primary_pars, env_interpolator, secondary_models) {
  
  ## Get the environmental conditions
  
  env_conditions <-  map(env_interpolator, ~.(time))
  
  ## Apply the secondary models
  
  p <- apply_secondary_models(env_conditions, secondary_models)
  
  state <- as.list(state)
  
  # dC <- -log(10)*state$C/p$D
  
  alpha <- 1
  beta <- (1 - p$Nres/state$N)
  
  dN <- -log(10)*state$N/p$D * alpha * beta
  
  res <- c(dN)
  return(list(res))
  
}

#' Geeraerd model without shoulder for predictions under dynamic conditions parameterized on k
#' 
#' The functions is defined to be called from [ode()] from the deSolve package.
#' 
#' @param time a 1D numeric with the treatment time
#' @param state a named numeric vector with two element defining the microbial concentration and 
#' the value of the ideal substance 'C'
#' @param primary_pars ignored. Just a dumpster for the param argument from ode
#' @param env_interpolator  an interpolator for the environmental conditions as returned by [approx_env()]
#' @param secondary_models a nested list defining the secondary models as per [apply_secondary_models()]
#' @param tol_time a tolerance for the time (to avoid singularities for t=0). See Details
#' 
#' @returns a list with two elements: dN/dt and dC/dt, as per [ode()]
#'
#'
dyna_Geeraerd_noSL_k <- function(time, state, primary_pars, env_interpolator, secondary_models) {
  
  ## Get the environmental conditions
  
  env_conditions <-  map(env_interpolator, ~.(time))
  
  ## Apply the secondary models
  
  p <- apply_secondary_models(env_conditions, secondary_models)
  
  state <- as.list(state)
  
  # dC <- - p$k * state$C
  
  alpha <- 1
  beta <- (1 - p$Nres/state$N)
  
  dN <- - p$k * state$N * alpha * beta
  
  res <- c(dN)
  return(list(res))
  
}








