
#'
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

#'
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








