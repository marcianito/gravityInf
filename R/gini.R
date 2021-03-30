#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2018), mreich@@posteo.de
#' @import 
#' @examples missing

gini = function(
    obs,
    sim
){
  ## debugging
  # obs = c(0,1,2,3,4,5)
  # sim = c(0,1,2,3,4,5)
  # for univariate time steps
  timestep = 1
  ## initiate lorenz curve results
  # for "complete" distribution
  Ages = 0
  # for observed distribution
  Alorenz = 0
  # calculate complete area
  for(i in 1:(length(obs) - 1)){
    # i = 1
    Ai = (obs[i+1] - obs[i] ) * timestep * 0.5 + obs[i] * timestep
    Ages = Ages + Ai
  }
  # calculate area below observations
  for(u in 1:(length(sim) - 1)){
    # u = 1
    Au = (sim[u+1] - sim[u] ) * timestep * 0.5 + sim[u] * timestep
    Alorenz = Alorenz + Au
  }
  # calculate gini coefficient
  gini_coef = (Ages - Alorenz) / Ages
  # return value
  return(gini_coef)
}

