# - ? 
NBapproxVARHSMM_generate <- function(N, parms, scale = TRUE, plt = TRUE) 
{
  
  K <- length(parms$Sigma)
  D <- nrow(parms$Sigma[[1]])
  P <- dim(parms$beta[[1]])[3]
  
  beta_0 <- parms$beta_0
  beta <- parms$beta
  Sigma <- parms$Sigma
  gamma <- parms$gamma
  delta <- parms$delta
  if(is.null(delta))  delta <- solve(t(diag(K)-gamma + 1), rep(1,K))
  lambda <- parms$lambda
  phi <- parms$phi
  
  
  state_seq <- numeric(N)
  obs <- matrix(NA, nrow = N, ncol = D)
  
  
  state_seq[1:P] <- sample(1:K, 1, prob = delta)
  obs[1:P, ] <- rmvnorm(P, beta_0[[state_seq[1]]], 
                        Sigma[[state_seq[1]]])
  
  state = state_seq[P]
  start = P + 1
  
  while(any(state_seq == 0)) {
    
    
    duration = 1 + rnbinom(1, mu = lambda[state], size = phi[state])
    end = start + duration - 1
    
    if (end > N) {
      end = N
      duration = length(start:end)
    } 
    
    state_seq[start:end] = rep(state, duration)
    
    for (n in start:end) {
      temp <- rep(0, D)
      for (p in 1:P) {
        temp = temp + (beta[[state]][, , p] %*% obs[n-p, ])
      }
      obs[n, ] <- rmvnorm(1, beta_0[[state]] + temp, Sigma[[state]])
    }
    start = end + 1
    probs = gamma[state_seq[end], ]
    state = sample(1:K, 1, prob = probs)
  }
  
  if (scale) {
    #obs <- standardize(obs)
    col_m <- colMeans(obs)
    col_s <- sqrt(colVars(obs))
    obs <- (obs - matrix(col_m, nrow = N, ncol = D, byrow = TRUE))/matrix(col_s, nrow = N, ncol = D, byrow = TRUE)
  } 
  
  
  if (plt) {
    
    # - multivariate single plot
    par(mfrow = c(1, 1))
    plot(1:N, type = "n",  ylab = "", xlab = "Time",
         cex = 0.5, las = 1,
         xlim = c(1, N),
         ylim = c(min(obs) - 1, max(obs) + 1), cex.lab = 1.1)
    if (K>2) {
      z_col <- brewer.pal(n = K, name = "Spectral")
    } else {
      z_col <- c("lightsalmon", "lightblue1")
    }
    
    z_aux <- sapply(1:N, function(n) z_col[state_seq[n]])
    cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
    col_all <- rainbow(D)
    for (j in 1:(length(cp_loc) - 1)) { 
      rect(cp_loc[j], ceiling(min(obs) - 2) , cp_loc[j+1], ceiling(max(obs) + 2) , 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
    }
    
    for (d in 1:D) {
      lines(obs[, d], pch = 20,
            col = col_all[d], cex = 1, type = "o")
    }
    
    # - separate univariate plots
    readline(prompt="Plot, univariate: Press [enter] to continue")
    n_row = 5
    n_col = ceiling(D/n_row)
    par(mfrow = c(n_row, n_col))
    
    par(mfrow = c(n_row, n_col), 
        mai = c(0.4, 0.55, 0.1, 0.3))
    
    for (dd in 1:D) {
      obs_univ <- obs[, dd]
      plot.new()
      plot.window(ylim = c(min(obs_univ) - 0.3, max(obs_univ) + 0.3), 
                  xlim = c(1, N),  xlab = "Time",  ylab = paste("y_", dd, sep=""), 
                  cex.lab = 1.1)
      title(xlab = "Time", ylab = paste("y_", dd, sep=""))
      for (j in 1:(length(cp_loc) - 1)) { 
        rect(cp_loc[j], ceiling(min(obs_univ)-2), cp_loc[j+1], ceiling(max(obs_univ) + 2), 
             col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
      }
      lines(1:N, obs_univ, pch = 16, cex = 1, col = col_all[dd], type = "o")
      axis(1); axis(2)
    }
  }
  
  if (scale) {
    return(list(obs = obs, state = state_seq, col_m = col_m, col_s = col_s))
  }
  else{
    return(list(obs = obs, state = state_seq))
  }
}

# - ? 
NBapproxVARHSMM_plot <- function(obs, state) 
{
  
  D <- ncol(obs)
  N <- nrow(obs)
  K <- length(unique(state))
  state_seq <- state
  
  # - multivariate single plot
  par(mfrow = c(1, 1))
  plot(1:N, type = "n",  ylab = "", xlab = "Time",
       cex = 0.5, las = 1,
       xlim = c(1, N),
       ylim = c(min(obs) - 1, max(obs) + 1), cex.lab = 1.1)
  if (K>2) {
    z_col <- brewer.pal(n = K, name = "Spectral")
  } else {
    z_col <- c("lightsalmon", "lightblue1")
  }
    
  z_aux <- sapply(1:N, function(n) z_col[state_seq[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
  col_all <- rainbow(D)
  for (j in 1:(length(cp_loc) - 1)) { 
    rect(cp_loc[j], ceiling(min(obs) - 2) , cp_loc[j+1], ceiling(max(obs) + 2) , 
         col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
  }
    
  for (d in 1:D) {
    lines(obs[, d], pch = 20,
          col = col_all[d], cex = 1, type = "o")
  }
    
  # - separate univariate plots
  #readline(prompt="Plot, univariate: Press [enter] to continue")
  n_row = 5
  n_col = ceiling(D/n_row)
  par(mfrow = c(n_row, n_col))
    
  par(mfrow = c(n_row, n_col), 
      mai = c(0.4, 0.55, 0.1, 0.3))
    
  for (dd in 1:D) {
    obs_univ <- obs[, dd]
    plot.new()
    plot.window(ylim = c(min(obs_univ) - 0.3, max(obs_univ) + 0.3), 
                xlim = c(1, N),  xlab = "Time",  ylab = paste("y_", dd, sep=""), 
                cex.lab = 1.1)
    title(xlab = "Time", ylab = paste("y_", dd, sep=""))
    for (j in 1:(length(cp_loc) - 1)) { 
      rect(cp_loc[j], ceiling(min(obs_univ)-2), cp_loc[j+1], ceiling(max(obs_univ) + 2), 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
    }
    lines(1:N, obs_univ, pch = 16, cex = 1, col = col_all[dd], type = "o")
    axis(1); axis(2)
  }
  
  return()
}

#  -
NBapproxVARHSMM_getPredictive <- function(fit, m, obs, pseudo = FALSE, L1_ball = TRUE, ndraw = 50)
{
  sims <- rstan::extract(fit)
  niter <- nrow(sims$beta)
  if(ndraw > niter) {stop("ndraw is larger than the number of  mcmc iterations:
                           choose a value between 1 and niter")}
  D <- ncol(obs)
  N <- nrow(obs)
  y_hat <- array(NA, dim = c(ndraw, N, D))
  log_pred_density <- array(NA, dim = c(ndraw, N))
  z_star <- array(NA, dim = c(ndraw, N))
  
  gamma_sims <- get_gamma_sims(sims$gamma)
  beta_0_sims <- sims$beta_0
  beta_sims <- sims$beta
  theta_sims <- get_theta_sims(sims, L1_ball)
  tau_sims <- sims$tau
  mu_sims <- get_mu_sims(sims, obs, L1_ball)
  Sigma_sims <- get_Sigma_sims(sims, pseudo)
  lambda_sims <- sims$lambda
  ## Cuase I goofed in the stan code 
  try(phi_sims <- 1/sims$phim1)
  if(length(phi_sims)==0){
    phi_sims <- 1/exp(sims$log_phim1)
  }
  K <- dim(gamma_sims)[2]
  
  
  beta_0_hat <- lapply(1:K, function(kk) { 
    apply(beta_0_sims[, kk, ], 2, mean)})
  beta_hat <- lapply(1:K, function(kk) { 
    apply(beta_sims[, kk, ], 2, mean)})
  theta_hat <- lapply(1:K, function(kk) { 
    apply(theta_sims[, kk, ], 2, mean)})
  tau_hat <- lapply(1:K, function(kk) { 
    apply(tau_sims[, kk, ], 2, mean)})
  gamma_hat <- apply(gamma_sims,  c(1, 2), mean)
  Sigma_hat <- lapply(1:K, function(kk) {
    apply(Sigma_sims[, kk, , ], c(2, 3), mean)})
  lambda_hat <- colMeans(lambda_sims)
  phi_hat <- colMeans(phi_sims)
  
  
  parms_hat <- list()
  parms_hat$beta_0 <- beta_0_hat
  parms_hat$beta <- beta_hat
  parms_hat$theta <- theta_hat
  parms_hat$Sigma <- Sigma_hat
  parms_hat$tau <- tau_hat
  parms_hat$gamma <- gamma_hat
  parms_hat$lambda <- lambda_hat
  parms_hat$phi <- phi_hat
  
  
  idxs <- sample(1:niter, ndraw, replace = FALSE)
  i <- 1
  cat(" ... sampling from posterior predictive... \n")
  
  for (t in idxs) {
    if((i %% (ndraw/10)) == 0) {
      cat(" ...", as.integer((i/ndraw)*100), "% \n")
    }
    
    beta_0_temp <- lapply(1:K, function(k, t) {beta_0_sims[t, k, ]}, t)
    beta_temp <-  lapply(1:K, function(k, t) {beta_sims[t, k, ]}, t)
    theta_temp <- lapply(1:K, function(k, t) {theta_sims[t, k, ]}, t)
    tau_temp <- lapply(1:K, function(k, t) {tau_sims[t, k, ]}, t)
    gamma_temp <- gamma_sims[, , t]
    Sigma_temp <- lapply(1:K, function(k, t) {Sigma_sims[t, k, , ]}, t)
    lambda_temp <- lambda_sims[t, ]
    phi_temp <- phi_sims[t, ]
    
    parms <- list()
    parms$beta_0 <- beta_0_temp
    parms$beta <- beta_temp
    parms$theta <- theta_temp
    parms$Sigma <- Sigma_temp
    parms$tau <- tau_temp
    parms$gamma <- gamma_temp
    parms$lambda <- lambda_temp
    parms$phi <- phi_temp
    
    z_star[i,] <- NBapproxVARHSMM_Viterbi(obs, m, parms, pseudo, L1_ball, draw = FALSE)
    
    y_hat[i, , ] <- t(sapply(1:N, function(n) {rmvnorm(1, mu_sims[t, z_star[i,n], n, ], 
                                                       Sigma_temp[[z_star[i,n]]])}))
    
    log_pred_density[i,] <- t(sapply(1:N, function(n) {dmvnorm(obs[n,], mu_sims[t, z_star[i,n], n, ], 
                                                               Sigma_temp[[z_star[i,n]]], log = TRUE)}))
    
    i <- i + 1
  }
  

  z_hat <- NBapproxVARHSMM_Viterbi(obs, m, parms_hat, pseudo, L1_ball, 
                                   draw = FALSE) 
  
  return(list(y_hat = y_hat, 
              z_hat = z_hat,
              z_star = z_star,
              log_pred_density = apply(log_pred_density, 2, logSumExp) - log(ndraw)
              ))
}

#  -
NBapproxVARHSMM_getPredictive_sims <- function(sims, m, obs, pseudo = FALSE, L1_ball = TRUE, ndraw = 50)
{
  #sims <- rstan::extract(fit)
  niter <- nrow(sims$beta)
  if(ndraw > niter) {stop("ndraw is larger than the number of  mcmc iterations:
                           choose a value between 1 and niter")}
  D <- ncol(obs)
  N <- nrow(obs)
  y_hat <- array(NA, dim = c(ndraw, N, D))
  log_pred_density <- array(NA, dim = c(ndraw, N))
  
  gamma_sims <- get_gamma_sims(sims$gamma)
  beta_0_sims <- sims$beta_0
  beta_sims <- sims$beta
  theta_sims <- get_theta_sims(sims, L1_ball)
  tau_sims <- sims$tau
  mu_sims <- get_mu_sims(sims, obs, L1_ball)
  Sigma_sims <- get_Sigma_sims(sims, pseudo)
  lambda_sims <- sims$lambda
  ## Cuase I goofed in the stan code 
  try(phi_sims <- 1/sims$phim1)
  if(length(phi_sims)==0){
    phi_sims <- 1/exp(sims$log_phim1)
  }
  K <- dim(gamma_sims)[2]
  
  
  beta_0_hat <- lapply(1:K, function(kk) { 
    apply(beta_0_sims[, kk, ], 2, mean)})
  beta_hat <- lapply(1:K, function(kk) { 
    apply(beta_sims[, kk, ], 2, mean)})
  theta_hat <- lapply(1:K, function(kk) { 
    apply(theta_sims[, kk, ], 2, mean)})
  tau_hat <- lapply(1:K, function(kk) { 
    apply(tau_sims[, kk, ], 2, mean)})
  gamma_hat <- apply(gamma_sims,  c(1, 2), mean)
  Sigma_hat <- lapply(1:K, function(kk) {
    apply(Sigma_sims[, kk, , ], c(2, 3), mean)})
  lambda_hat <- colMeans(lambda_sims)
  phi_hat <- colMeans(phi_sims)
  
  
  parms_hat <- list()
  parms_hat$beta_0 <- beta_0_hat
  parms_hat$beta <- beta_hat
  parms_hat$theta <- theta_hat
  parms_hat$Sigma <- Sigma_hat
  parms_hat$tau <- tau_hat
  parms_hat$gamma <- gamma_hat
  parms_hat$lambda <- lambda_hat
  parms_hat$phi <- phi_hat
  
  
  idxs <- sample(1:niter, ndraw, replace = FALSE)
  i <- 1
  cat(" ... sampling from posterior predictive... \n")
  
  for (t in idxs) {
    if((i %% (ndraw/10)) == 0) {
      cat(" ...", as.integer((i/ndraw)*100), "% \n")
    }
    
    beta_0_temp <- lapply(1:K, function(k, t) {beta_0_sims[t, k, ]}, t)
    beta_temp <-  lapply(1:K, function(k, t) {beta_sims[t, k, ]}, t)
    theta_temp <- lapply(1:K, function(k, t) {theta_sims[t, k, ]}, t)
    tau_temp <- lapply(1:K, function(k, t) {tau_sims[t, k, ]}, t)
    gamma_temp <- gamma_sims[, , t]
    Sigma_temp <- lapply(1:K, function(k, t) {Sigma_sims[t, k, , ]}, t)
    lambda_temp <- lambda_sims[t, ]
    phi_temp <- phi_sims[t, ]
    
    parms <- list()
    parms$beta_0 <- beta_0_temp
    parms$beta <- beta_temp
    parms$theta <- theta_temp
    parms$Sigma <- Sigma_temp
    parms$tau <- tau_temp
    parms$gamma <- gamma_temp
    parms$lambda <- lambda_temp
    parms$phi <- phi_temp
    
    z_star <- NBapproxVARHSMM_Viterbi(obs, m, parms, pseudo, L1_ball, draw = FALSE)
    
    y_hat[i, , ] <- t(sapply(1:N, function(n) {rmvnorm(1, mu_sims[t, z_star[n], n, ], 
                                                       Sigma_temp[[z_star[n]]])}))
    
    log_pred_density[i,] <- t(sapply(1:N, function(n) {dmvnorm(obs[n,], mu_sims[t, z_star[n], n, ], 
                                                               Sigma_temp[[z_star[n]]], log = TRUE)}))
    
    i <- i + 1
  }
  
  
  z_hat <- NBapproxVARHSMM_Viterbi(obs, m, parms_hat, pseudo, L1_ball, 
                                   draw = FALSE) 
  
  return(list(y_hat = y_hat, 
              z_hat = z_hat,
              log_pred_density = apply(log_pred_density, 2, logSumExp) - log(ndraw)
              ))
}

## This function does a one-step ahead prediction ahead prediction 
# conditioning on the previous observation and current state according to the viterbi
# most likely - but this is integreated over the posterior
## The full Bayesian predictive would rather than taking the most likely state, 
# would integrate over the dist of the states (as well as the parameters - cause he samples from the posterior)
# so all we would need to change would be to calculate the previous states probabilities 
# also integrate over state uncertainty - i would run the bayesian get state probs first i guess, 
# would need to make sure it didn't use the future - i think at the moment using the viterbi cheats slightly
## Note the z_hat's are just from a posterior mean though!

## First I will edit the original function to calculate the predictive density



# - 
NBapproxVARHSMM_Viterbi <- function(obs, m, parms, pseudo = FALSE, 
                                    L1_ball = TRUE, draw = FALSE)
{
  
  N <- nrow(obs)
  D <- ncol(obs)
  K <- length(m)
  M <- sum(m)
  
  gamma <- parms$gamma
  lambda <- parms$lambda
  phi <- parms$phi
  
  B <- B_matrix_negativeBinomial(K, m_vect = m, a_mat = gamma, lambda, phi)
  allprobs_aggr <- VAR_emissions_aggr(obs, m, parms, pseudo, L1_ball)
  xi <- matrix(0, N, M)
  delta <- solve(t(diag(M)-B + 1), rep(1,M))
  foo <- delta * allprobs_aggr[1, ]
  xi[1, ] <- foo/sum(foo)
  
  for (t in 2:N) {
    foo <- apply(xi[t-1, ] * B, 2, max) * allprobs_aggr[t, ]
    xi[t, ] <- foo/sum(foo)
  }
  z_star <- numeric(N)
  
  # - Drawing 
  if (draw == TRUE) {
    z_star[N] <- sample(1:M, size = 1, prob = xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- sample(1:M, size = 1, prob = B[, z_star[t+1]] * xi[t, ])
    }
  } else {
    # - Maximizing
    z_star[N] <- which.max(xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- which.max(B[, z_star[t+1]] * xi[t, ])
    }
  }
  
  # standard state sequence (not aggregate)
  cm = c(0, cumsum(m))
  for (j in 2:(K+1)) {
    idx = (cm[j - 1]+1):cm[j]
    z_star[z_star %in% idx] = j - 1
  }
  
  return(z_star)
}


# - 
NBapproxVARHSMM_getBayesEstimates <- function(sims, pseudo = FALSE, L1_ball = TRUE)
{
  gamma_sims <- get_gamma_sims(sims$gamma)
  beta_0_sims <- sims$beta_0
  beta_sims <- sims$beta
  theta_sims <- get_theta_sims(sims, L1_ball)
  
  
  tau_sims <- sims$tau
  Sigma_sims <- get_Sigma_sims(sims, pseudo)
  lambda_sims <- sims$lambda
  ## Cause I goofed in the stan code 
  try(phi_sims <- 1/sims$phim1)
  if(length(phi_sims)==0){
    phi_sims <- 1/exp(sims$log_phim1)
  }
  K <- dim(gamma_sims)[2]
  
  beta_0_hat <- lapply(1:K, function(kk) { 
    apply(beta_0_sims[, kk, ], 2, mean)})
  beta_hat <- lapply(1:K, function(kk) { 
    apply(beta_sims[, kk, ], 2, mean)})
  tau_hat <- lapply(1:K, function(kk) { 
    apply(tau_sims[, kk, ], 2, mean)})
  gamma_hat <- apply(gamma_sims,  c(1, 2), mean)
  Sigma_hat <- lapply(1:K, function(kk) {
    apply(Sigma_sims[, kk, , ], c(2, 3), mean)})
  lambda_hat <- colMeans(lambda_sims)
  phi_hat <- colMeans(phi_sims)
  
  
  parms_hat <- list()
  parms_hat$beta_0 <- beta_0_hat
  parms_hat$beta <- beta_hat
  parms_hat$Sigma <- Sigma_hat
  parms_hat$tau <- tau_hat
  parms_hat$gamma <- gamma_hat
  parms_hat$lambda <- lambda_hat
  parms_hat$phi <- phi_hat
  
  if (L1_ball) {
    theta_hat <- lapply(1:K, function(kk) { 
      apply(theta_sims[, kk, ], 2, mean)})
    parms_hat$theta <- theta_hat
  }
  
  return(parms_hat)
}

# - 
NBapproxVARHSMM_forwardBackwards <- function(obs, m, allprobs_aggr, parms)
{
  
  #lamdab <- parms$lambda
  lambda <- parms$lambda
  phi <- parms$phi
  gamma <- parms$gamma

  
  M <- sum(m)
  N <- nrow(obs)
  lalpha <- lbeta <- matrix(NA, M, N)
  
  B <- B_matrix_negativeBinomial(K, m_vect = m, a_mat = gamma, lambda, phi)
  delta <- rep(1, M)
  
  foo <- delta*allprobs_aggr[1, ]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[, 1] <- log(foo) + lscale
  
  
  for (i in 2:N) 
  {
    foo <- foo%*%B*allprobs_aggr[i,]
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
    lalpha[, i] <- log(foo) +lscale
  }
  llk <- lscale
  
  lbeta[, N] <- rep(0, M)
  foo <- rep(1/M, M)
  lscale <- log(M)
  
  for (i in (N-1):1) {
    foo <- B%*%(allprobs_aggr[i+1,]*foo)
    lbeta[,i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  return(list(la = lalpha, lb = lbeta, llk = llk))
}

# - 
NBapproxVARHSMM_stateProbs <- function(sims, obs, m, pseudo = FALSE, L1_ball = TRUE)
{
  parms_hat <- NBapproxVARHSMM_getBayesEstimates(sims, pseudo, L1_ball)
  M <- sum(m)
  N <- nrow(obs)
  
  allprobs_aggr <- VAR_emissions_aggr(obs, m, parms_hat, pseudo, L1_ball)
  fb <- NBapproxVARHSMM_forwardBackwards(obs, m, allprobs_aggr, parms_hat)
  la <- fb$la
  lb <- fb$lb
  llk <- fb$llk
  
  stateprobs <- matrix(NA, ncol = N, nrow = M)
  for (n in 1:N) stateprobs[, n] <- exp(la[, n] + lb[, n] - llk)
  return(stateprobs)
}

# The previous version of this function just does so for the posterior means, 
# why not use the full posterior sample 
# - 
NBapproxVARHSMM_stateProbs_Bayesian <- function(sims, obs, m, pseudo = FALSE, L1_ball = TRUE, ndraw = 50)
{
  #sims <- rstan::extract(fit)
  niter <- nrow(sims$beta)
  if(ndraw > niter) {stop("ndraw is larger than the number of  mcmc iterations:
                           choose a value between 1 and niter")}
  
  gamma_sims <- get_gamma_sims(sims$gamma)
  beta_0_sims <- sims$beta_0
  beta_sims <- sims$beta
  tau_sims <- sims$tau
  Sigma_sims <- get_Sigma_sims(sims, pseudo)
  theta_sims <- get_theta_sims(sims, L1_ball)
  lambda_sims <- sims$lambda
  phi_sims <- 1/sims$phim1
  K <- dim(gamma_sims)[2]
  N <- nrow(obs)
  M <- sum(m)
  
  N_MCMC <- length(sims$lp__)
  stateprobs <- array(NA, dim = c(ndraw, M, N))

  #parms_hat <- VARHMM_getBayesEstimates(sims, pseudo, L1_ball)
  #for(i in 1:N_MCMC){
  idxs <- sample(1:niter, ndraw, replace = FALSE)
  i <- 1
  cat(" ... sampling from posterior predictive... \n")
  
  for (t in idxs) {
    if((i %% (ndraw/10)) == 0) {
      cat(" ...", as.integer((i/ndraw)*100), "% \n")
    }
    parms_hat <- list()
    parms_hat$beta_0 <- lapply(1:K, function(kk) {beta_0_sims[t,kk,]})
    parms_hat$beta <- lapply(1:K, function(kk) {beta_sims[t,kk,]})
    parms_hat$Sigma <- lapply(1:K, function(kk) {Sigma_sims[t,kk,,]})
    parms_hat$tau <- lapply(1:K, function(kk) {tau_sims[t,kk,]})
    parms_hat$gamma <- gamma_sims[,,t]
    parms_hat$lambda <- lambda_sims[t,]
    parms_hat$phi <- phi_sims[t,]
    if (L1_ball) {
      parms_hat$theta <- lapply(1:K, function(kk) {theta_sims[t,kk,]}) 
    }
    
    allprobs_aggr <- VAR_emissions_aggr(obs, m, parms_hat, pseudo, L1_ball)
    fb <- NBapproxVARHSMM_forwardBackwards(obs, m, allprobs_aggr, parms_hat)
    la <- fb$la
    lb <- fb$lb
    llk <- fb$llk
    
    #stateprobs <- matrix(NA, ncol = N, nrow = K)
    for (n in 1:N) stateprobs[i,, n] <- exp(la[, n] + lb[, n] - llk)
    
    i <- i + 1
    
  }
  
  
  return(apply(stateprobs, c(2, 3), mean))
}

# - 
# - time varying correlation
NBapproxVARHSMM_getCorrelation <- function(fit, obs, m, dimensions, 
                                           n_draw = 200, 
                                           pseudo, L1_ball )
{
  sims <- rstan::extract(fit)
  gamma_sims <- get_gamma_sims(sims$gamma)
  beta_0_sims <- sims$beta_0
  beta_sims <- sims$beta
  theta_sims <- get_theta_sims(sims, L1_ball)
  tau_sims <- sims$tau
  mu_sims <- get_mu_sims(sims, obs, L1_ball)
  Sigma_sims <- get_Sigma_sims(sims, pseudo)
  lambda_sims <- sims$lambda
  phi_sims <- 1/sims$phim1
  
  K <- dim(gamma_sims)[2]
  n_MCMC <- dim(gamma_sims)[3]
  D <- ncol(obs)
  N <- nrow(obs)
  
  corr_seq_sims <- matrix(NA, nrow=n_draw, ncol=N)
  corr_sims <- list()
  
  for (kk in 1:K) {
    temp <- numeric(n_MCMC)
    for (tt in 1:n_MCMC) {
      S <- Sigma_sims[tt, kk, , ]
      R <- cov2cor(S)
      temp[tt] <- R[dimensions[1], dimensions[2]]
    }
    corr_sims[[kk]] <- temp
  }
  
  idxs <- sample(1:n_MCMC, n_draw, replace = FALSE)
  i <- 1
  cat(" ... predictive time-varying correlation.. \n")
  
  for (t in idxs) {
    if((i %% (n_draw/10)) == 0) {
      cat(" ...", as.integer((i/n_draw)*100), "% \n")
    }
    
    beta_0_temp <- lapply(1:K, function(k, t) {beta_0_sims[t, k, ]}, t)
    beta_temp <-  lapply(1:K, function(k, t) {beta_sims[t, k, ]}, t)
    theta_temp <- lapply(1:K, function(k, t) {theta_sims[t, k, ]}, t)
    tau_temp <- lapply(1:K, function(k, t) {tau_sims[t, k, ]}, t)
    gamma_temp <- gamma_sims[, , t]
    Sigma_temp <- lapply(1:K, function(k, t) {Sigma_sims[t, k, , ]}, t)
    lambda_temp <- lambda_sims[t, ]
    phi_temp <- phi_sims[t, ]
    
    parms <- list()
    parms$beta_0 <- beta_0_temp
    parms$beta <- beta_temp
    parms$theta <- theta_temp
    parms$Sigma <- Sigma_temp
    parms$tau <- tau_temp
    parms$gamma <- gamma_temp
    parms$lambda <- lambda_temp
    parms$phi <- phi_temp
    
    z_hat <- NBapproxVARHSMM_Viterbi(obs, m, parms, pseudo, L1_ball, draw = FALSE)
    
    for (n in 1:N){
      k <- z_hat[n]
      corr_seq_sims[i, n] <- corr_sims[[k]][t]
    }
    i <- i + 1
  }
  
  return(corr_seq_sims)
}


# - 
NBapproxVARHSMM_plotStateProbs <- function(fit, obs, m, pseudo = FALSE, L1_ball = TRUE, switch = FALSE)
{
  sims <- rstan::extract(fit)
  K <- dim(sims$gamma)[2]
  #M <- length(m)
  M <- sum(m)
  N <- nrow(obs)
  time <- 1:N
  if (K>2) {
    z_col <- brewer.pal(n = K, name = "Spectral")
  } else {
    if (switch) {
      z_col <- c("lightblue1", "lightsalmon")
      
    } else {
      z_col <- c("lightsalmon", "lightblue1")
      
    }
  }
  
  state_probs <- NBapproxVARHSMM_stateProbs(sims, obs, m, pseudo, L1_ball)
  parms_hat <- NBapproxVARHSMM_getBayesEstimates(sims, pseudo, L1_ball)
  z_hat <- NBapproxVARHSMM_Viterbi(obs, m, parms_hat, pseudo, L1_ball, draw = FALSE) 
  z_aux <- sapply(1:N, function(n) z_col[z_hat[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
  
  # retrieving state probs from augumentated matrix
  state_probs_mod <- matrix(NA, K, N)
  for (nn in 1:N) {
    for (kk in 1:K) {
      if (kk == 1) {
        start = 1
        end = m[kk]
      } else {
        start = sum(m[1:(kk-1)]) + 1
        end = sum(m[1:kk])
      }
      state_probs_mod[kk, nn] = sum(state_probs[start:end, nn])
    }
  }
  
  
  plot(1:N, type = "n", xlab = "Time", ylab = "Prob state", 
       cex = 0.5, las = 1,
       xlim = c(1, N),
       ylim = c(0,1), cex.lab = 1.1)
  
  plot_p <- matrix(NA,nrow=N-1,ncol=2*K)
  a <- 0
  
  for (n in 2:N) {
    
    for (j in 1:K) {
      plot_p[n-1, (j*2-1)] <- state_probs_mod[j, n-1]
      plot_p[n-1, j*2] <- state_probs_mod[j, n]
      if	(j==1){col_states<- z_col[1]}					
      if	(j==2){col_states<- z_col[2]}					
      if	(j==3){col_states<- z_col[3]}	
      if	(j==4){col_states<- z_col[4]}	
      
      # (at some point need to make function to 
      # generalize for all K, like this is redundant)
      if	(j==1){				
        point_1<-a
        point_2<-point_1+plot_p[n-1,(j*2-1)]
        point_4<-a
        point_3<-point_4+plot_p[n-1,(j*2)]	}
      
      if	(j==2){				
        point_1<-a+plot_p[n-1,(j-1)*2-1]
        point_2<-point_1+plot_p[n-1,(j*2-1)]
        point_4<-a+plot_p[n-1,(j-1)*2]
        point_3<-point_4+plot_p[n-1,(j*2)]	}
      
      if	(j==3){				
        point_1<-a+plot_p[n-1,(j-2)*2-1]+plot_p[n-1,(j-1)*2-1]
        point_2<-point_1+plot_p[n-1,(j*2-1)]
        point_4<-a+plot_p[n-1,(j-2)*2]+plot_p[n-1,(j-1)*2]
        point_3<-point_4+plot_p[n-1,(j*2)]}
      if (j==4) {
        point_1 <- a+ plot_p[n-1,(j-3)*2-1] + plot_p[n-1,(j-2)*2-1]+plot_p[n-1,(j-1)*2-1]
        point_2 <- point_1+plot_p[n-1,(j*2-1)]
        point_4 <- a + plot_p[n-1,(j-3)*2] + plot_p[n-1,(j-2)*2]+plot_p[n-1,(j-1)*2]
        point_3 <- point_4+plot_p[n-1,(j*2)]}
      
      polygon(c(time[n-1],time[n-1],time[n],time[n]),
              c(point_1,point_2,point_3,point_4),col=scales::alpha(col_states, 0.8), border=NA)
      lines(c(time[n-1],time[n]),c(point_2,point_3), col=scales::alpha(col_states, 0.5))
    }
  }
}


## Hazard function for the Negative-Binomail
c_hazard_NegBinom <- function(r, lambda, phi){
  ## Shifted so that the minimum value is 1
  return(dnbinom(r-1, mu = lambda, size = phi)/(1-pnbinom(r-2, mu = lambda, size = phi)))
}


# B_ii_negativeBinomial;
B_ii_negativeBinomial <- function(m_i, lambda_i, phi_i){
  B_ii <- matrix(0, nrow = m_i, ncol = m_i);
  if(m_i > 1){
    for(i in 1:(m_i-1)){
      B_ii[i, i+1] = 1 - c_hazard_NegBinom(i, lambda_i, phi_i);
    }
  }
  B_ii[m_i, m_i] = 1 - c_hazard_NegBinom(m_i, lambda_i, phi_i);
  return(B_ii)
}

# B_ij_negativeBinomial;
B_ij_negativeBinomial <- function(m_i, m_j, a_ij, lambda_i, phi_i){
  B_ij <- matrix(0, nrow = m_i, ncol = m_j);
  for( i in 1:m_i){
    B_ij[i, 1] =  a_ij * c_hazard_NegBinom(i, lambda_i, phi_i);
  }
  return(B_ij)
}

# B_matrix_negativeBinomial;
B_matrix_negativeBinomial <- function(K, m_vect, a_mat, lambda_vec, phi_vec){
  sum_m = sum(m_vect);
  m_vect_temp <- rep(NA, K+1)
  B = matrix(0, nrow = sum_m, ncol = sum_m);
  m_vect_temp[1] = 0;
  for(i in 1:K){
    m_vect_temp[i+1] = m_vect[i];
  }
  for(i in 1:K){
    for(j in 1:K){
      if(i ==j){
        B[(sum(m_vect_temp[1:i])+1):sum(m_vect_temp[1:(i+1)]), (sum(m_vect_temp[1:j])+1):sum(m_vect_temp[1:(j+1)])] =
          B_ii_negativeBinomial(m_vect_temp[i+1], lambda_vec[i], phi_vec[i]);
      } 
      else{
        B[(sum(m_vect_temp[1:i])+1):sum(m_vect_temp[1:(i+1)]), (sum(m_vect_temp[1:j])+1):sum(m_vect_temp[1:(j+1)])] =
          B_ij_negativeBinomial(m_vect_temp[i+1], m_vect_temp[j+1], a_mat[i,j], lambda_vec[i], phi_vec[i]);
      }
    }
  }
  return(B)
}