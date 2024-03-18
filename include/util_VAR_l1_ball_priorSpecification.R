### VAR PARAMETERS ###

theta_prior_gen_l1_ball <- function(p, r_alpha, beta_gen, N_MC = 5000){
  # r ~ Exp(r_alpha)
  theta <- matrix(NA, nrow = N_MC, ncol = p)
  for(j in 1:N_MC){
    r <- rexp(1, r_alpha)
    beta <- beta_gen(p)
    theta[j, ] <- l1_ball_projection(beta, r)
  }
  return(theta)
}

## Requires:
## util.R/get_companion_matrix written by beniamo
## Simialr to:
## util.R/generate_sparse_beta

## mimicing l1-ball but with given sparsity level 
sparse_VAR_DExp_stability_prob <- function(N_MC, D, P, sparsity_level, lambda) 
{
  # N_MC is th enumber of monte-carlo sample sused to estimate the probability
  # D is the dimension of the data
  # P is the number AR lag
  # sparsity level is the targetted sparsity level (could potentiall extend to different sparisy tin different lags etc.)
  # lambda is the scale of the theta's
  
  theta <- array(NA, dim = c(D, D, P))
  #zeros <- sample(c(0, 1), P*D^2, replace = TRUE,
  #                prob = c(1 - sparsity_level, sparsity_level))
  stab <- array(NA, dim = c(N_MC))
  for(j in 1:N_MC){
    for (p in 1:P) {
      zeros <- sample(c(0, 1), D^2, replace = TRUE,
                      prob = c(1 - sparsity_level, sparsity_level))
      elements <- (1 - zeros) * rlaplace(D^2, 0, lambda)
      theta[, , p] <-  matrix(elements, D, D, byrow = T)
    }
    A_comp <- get_companion_matrix(theta)
    stab[j] <- all(abs(eigen(A_comp)$values) < 1)
  }
  
  return(stab)
}

sparse_VAR_Norm_stability_prob <- function(N_MC, D, P, sparsity_level, lambda) 
{
  # N_MC is th enumber of monte-carlo sample sused to estimate the probability
  # D is the dimension of the data
  # P is the number AR lag
  # sparsity level is the targetted sparsity level (could potentiall extend to different sparisy tin different lags etc.)
  # lambda is the scale of the theta's
  
  theta <- array(NA, dim = c(D, D, P))
  #zeros <- sample(c(0, 1), P*D^2, replace = TRUE,
  #                prob = c(1 - sparsity_level, sparsity_level))
  stab <- array(NA, dim = c(N_MC))
  for(j in 1:N_MC){
    for (p in 1:P) {
      zeros <- sample(c(0, 1), D^2, replace = TRUE,
                      prob = c(1 - sparsity_level, sparsity_level))
      elements <- (1 - zeros) * rnorm(D^2, 0, lambda)
      theta[, , p] <-  matrix(elements, D, D, byrow = T)
    }
    A_comp <- get_companion_matrix(theta)
    stab[j] <- all(abs(eigen(A_comp)$values) < 1)
  }
  
  return(stab)
}

## mimicing l1-ball but with given sparsity level 
non_sparse_VAR_DExp_stability_prob <- function(N_MC, D, P, lambda) 
{
  # N_MC is th enumber of monte-carlo sample sused to estimate the probability
  # D is the dimension of the data
  # P is the number AR lag
  # lambda is the scale of the theta's
  
  theta <- array(NA, dim = c(D, D, P))
  #zeros <- sample(c(0, 1), P*D^2, replace = TRUE,
  #                prob = c(1 - sparsity_level, sparsity_level))
  stab <- array(NA, dim = c(N_MC))
  for(j in 1:N_MC){
    for (p in 1:P) {
      elements <- rlaplace(D^2, 0, lambda)
      theta[, , p] <-  matrix(elements, D, D, byrow = T)
    }
    A_comp <- get_companion_matrix(theta)
    stab[j] <- all(abs(eigen(A_comp)$values) < 1)
  }
  
  return(stab)
}

## mimicing l1-ball but with given sparsity level 
non_sparse_VAR_Norm_stability_prob <- function(N_MC, D, P, lambda) 
{
  # N_MC is th enumber of monte-carlo sample sused to estimate the probability
  # D is the dimension of the data
  # P is the number AR lag
  # lambda is the scale of the theta's
  
  theta <- array(NA, dim = c(D, D, P))
  #zeros <- sample(c(0, 1), P*D^2, replace = TRUE,
  #                prob = c(1 - sparsity_level, sparsity_level))
  stab <- array(NA, dim = c(N_MC))
  for(j in 1:N_MC){
    for (p in 1:P) {
      elements <- rnorm(D^2, 0, lambda)
      theta[, , p] <-  matrix(elements, D, D, byrow = T)
    }
    A_comp <- get_companion_matrix(theta)
    stab[j] <- all(abs(eigen(A_comp)$values) < 1)
  }
  
  return(stab)
}


### DWELL DISTRIBUTIONS ###
# taken from Hadj-Amar, Jewson and Ficas (2022)

gamma_simplex_to_gamma_matrix <- function(gamma_simplex){
  p <- nrow(gamma_simplex)
  gamma_matrix <- matrix(0, nrow = p, ncol = p)
  for(i in 1:p){
    ind <- 0
    for(j in 1:p){
      if(i ==j){next}
      ind <- ind + 1
      gamma_matrix[i, j] <- gamma_simplex[i, ind]
    }
  }
  return(gamma_matrix)
}

geomDwellMean_mean_variance <- function(alpha_i, beta_i){
  return(list("mean" = (alpha_i + beta_i -1)/(beta_i -1), "var" = (alpha_i + beta_i -1)*(alpha_i + beta_i -2)/((beta_i -1)*(beta_i -2)) - ((alpha_i + beta_i -1)/(beta_i -1))^2))
}

geomDwell_mean_var_target_error <- function(alpha_i, beta_i, mean_target, var_target){
  temp <- geomDwellMean_mean_variance(alpha_i, beta_i)
  return(abs(mean_target - temp$mean) + abs(var_target - temp$var))
}

NBDwellMean_mean_variance <- function(a_0, b_0){## plus 1 to shift the support of the NB to >0
  return(list("mean" = a_0/b_0 + 1, "var" = a_0/(b_0^2)))
}

NBDwell_mean_var_target_error <- function(a_0, b_0, mean_target, var_target){
  temp <- NBDwellMean_mean_variance(a_0, b_0)
  return(abs(mean_target - temp$mean) + abs(var_target - temp$var))
}

NBDwellMean_IG_mean_variance <- function(a_0, b_0){## plus 1 to shift the support of the NB to >0
  return(list("mean" = b_0/(a_0 - 1) + 1, "var" = (b_0^2)/((a_0 - 1)^2*(a_0 - 2))))
}

NBDwell_IG_mean_var_target_error <- function(a_0, b_0, mean_target, var_target){
  temp <- NBDwellMean_IG_mean_variance(a_0, b_0)
  return(abs(mean_target - temp$mean) + abs(var_target - temp$var))
}


### Non-local Priors ###

dpMOM <- function(theta, tau, sigma2){
  return(theta^2/(tau*sigma2)*dnorm(theta, 0, sqrt(tau*sigma2)))
}

dpiMOM <- function(theta, tau, sigma2){
  return(sqrt(tau*sigma2)/(sqrt(pi)*theta^2)*exp(- tau*sigma2/theta^2))
}

dpeMOM <- function(theta, tau, sigma2){
  return(exp(sqrt(2) - (tau*sigma2)/(theta^2))*dnorm(theta, 0, sqrt(tau*sigma2)))
}

dpMOM_phi <- function(phi, tau, sigma2){
  return((-log(phi))^2/(tau*sigma2*phi)*dnorm(-log(phi), 0, sqrt(tau*sigma2)))
}

dpeMOM_phi <- function(phi, tau, sigma2){
  return(exp(sqrt(2) - (tau*sigma2)/((log(phi))^2))/phi*dnorm(log(phi), 0, sqrt(tau*sigma2)))
}

TVD_NegBinom_Geom <- function(p, rho, x_thresh){
  x_seq <- 1:x_thresh
  ## prob is the success probability 
  return(0.5*sum(abs(dgeom(x_seq - 1, prob = 1 - p) - dnbinom(x_seq - 1, mu = p/(1-p), size = rho))))
}

TVD_NegBinom_Geom_rhom1 <- function(p, rhom1, x_thresh){
  x_seq <- 1:x_thresh
  ## prob is the success probability 
  return(0.5*sum(abs(dgeom(x_seq - 1, prob = 1 - p) - dnbinom(x_seq - 1, mu = p/(1-p), size = 1/rhom1))))
}
