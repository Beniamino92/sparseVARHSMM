# - ? 
VAR_emissions_aggr <- function(obs, m, parms, pseudo = FALSE, L1_ball = TRUE)
{
  N <-dim(obs)[1]
  K <- length(m)
  
  # - emission probs
  allprobs <- VAR_emissions(obs, parms, pseudo, L1_ball)  
  
  # - emission aggregates
  allprobs_aggr = c()
  for (i in 1:K) {
    temp = matrix(rep(allprobs[, i], times = m[i]), nrow = N, ncol = m[i])
    allprobs_aggr = cbind(allprobs_aggr, temp)
  }
  allprobs_aggr
}


