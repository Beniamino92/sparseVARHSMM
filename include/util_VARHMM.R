# - 
VARHMM_generate <- function(N, parms, scale = FALSE, plt = TRUE) 
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
  
  # checking VAR stability, for each state
  for (kk in 1:K) {
    A_comp <- get_companion_matrix(beta[[kk]])
    if(!all(abs(eigen(A_comp)$values) < 1)) {
      cat("VAR[[", kk, "]]", " is not stable", "\n")
      stop()
    }
  }
  
  
  state_seq <- numeric(N)
  obs <- matrix(NA, nrow = N, ncol = D)
  state_seq[1:P] <- sample(1:K, 1, prob = delta)
  
  obs[1:P, ] <- rmvnorm(P, beta_0[[state_seq[1]]], 
                        Sigma[[state_seq[1]]])
  
  for (n in (P+1):N) {
    state <- sample(1:K, 1, prob = gamma[state_seq[n-1], ])
    state_seq[n] <- state
    temp <- rep(0, D)
    for (p in 1:P) {
      temp = temp + (beta[[state]][, , p] %*% obs[n-p, ])
    }
    obs[n, ] <- rmvnorm(1, beta_0[[state]] + temp, Sigma[[state]])
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
      rect(cp_loc[j], ceiling(min(obs)) - 1, cp_loc[j+1], ceiling(max(obs)) + 1, 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
    }
    
    for (d in 1:D) {
      lines(obs[, d], pch = 20,
            col = col_all[d], cex = 1, type = "o")
    }
    
    # - separate univariate plots
    invisible(readline(prompt="Plot, univariate: Press [enter] to continue"))
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
        rect(cp_loc[j], ceiling(min(obs)), cp_loc[j+1], ceiling(max(obs)), 
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


# - 
VARHMM_plot <- function(obs, state) 
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
    rect(cp_loc[j], ceiling(min(obs)) - 1, cp_loc[j+1], ceiling(max(obs)) + 1, 
         col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
  }
    
  for (d in 1:D) {
    lines(obs[, d], pch = 20,
          col = col_all[d], cex = 1, type = "o")
  }
    
  # - separate univariate plots
  #invisible(readline(prompt="Plot, univariate: Press [enter] to continue"))
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
      rect(cp_loc[j], ceiling(min(obs)), cp_loc[j+1], ceiling(max(obs)), 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
    }
    lines(1:N, obs_univ, pch = 16, cex = 1, col = col_all[dd], type = "o")
      axis(1); axis(2)
    }
  return()
}

# - 
VARHMM_viterbi <- function(obs, parms, pseudo = FALSE, L1_ball = TRUE, 
                            draw = FALSE) 
{
  beta_0 <- parms$beta_0
  beta <- parms$beta
  theta <- parms$theta
  tau <- parms$tau
  gamma <- parms$gamma  
  
  N <- nrow(obs)
  D <- ncol(obs)
  K <- length(beta_0)
  P_hat <- length(beta[[1]])
  P <- P_hat/(D*D)
  
  allprobs <- VAR_emissions(obs, parms, pseudo = pseudo, L1_ball = L1_ball)
  allprobs <- ifelse(!is.na(allprobs), allprobs, 1)
  
  gamma_init <- solve(t(diag(K)-gamma+1),rep(1,K))
  xi <- matrix(0, N, K)
  foo <- gamma_init * allprobs[1, ]
  xi[1, ] <- foo/sum(foo)
  
  for (t in 2:N) {
    foo <- apply(xi[t-1, ] * gamma, 2, max) * allprobs[t, ]
    xi[t, ] <- foo/sum(foo)
  }
  z_star <- numeric(N)
  
  # - Resampling
  if (draw == TRUE) {
    z_star[N] <- sample(1:K, size = 1, prob = xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- sample(1:K, size = 1, prob = gamma[, z_star[t+1]] * xi[t, ])
    }
  } else { # Maximizing
    z_star[N] <- which.max(xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- which.max(gamma[, z_star[t+1]] * xi[t, ])
    }
  }
  
  return(z_star)
}

# - 
VARHMM_get_predictive <- function(fit, obs, pseudo = FALSE, L1_ball = TRUE, 
                                   ndraw = 50)
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
  
  gamma_sims <- sims$gamma
  beta_0_sims <- sims$beta_0
  beta_sims <- sims$beta
  theta_sims <- get_theta_sims(sims, L1_ball)
  tau_sims <- sims$tau
  mu_sims <- get_mu_sims(sims, obs, L1_ball)
  Sigma_sims <- get_Sigma_sims(sims, pseudo)
  K <- dim(gamma_sims)[2]
  
  
  beta_0_hat <- lapply(1:K, function(kk) { 
    apply(beta_0_sims[, kk, ], 2, mean)})
  beta_hat <- lapply(1:K, function(kk) { 
    apply(beta_sims[, kk, ], 2, mean)})
  theta_hat <- lapply(1:K, function(kk) { 
    apply(theta_sims[, kk, ], 2, mean)})
  tau_hat <- lapply(1:K, function(kk) { 
    apply(tau_sims[, kk, ], 2, mean)})
  gamma_hat <- apply(gamma_sims,  c(2, 3), mean)
  Sigma_hat <- lapply(1:K, function(kk) {
    apply(Sigma_sims[, kk, , ], c(2, 3), mean)})
  
  
  parms_hat <- list()
  parms_hat$beta_0 <- beta_0_hat
  parms_hat$beta <- beta_hat
  parms_hat$theta <- theta_hat
  parms_hat$Sigma <- Sigma_hat
  parms_hat$tau <- tau_hat
  parms_hat$gamma <- gamma_hat
  
  
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
    gamma_temp <- gamma_sims[t, , ]
    Sigma_temp <- lapply(1:K, function(k, t) {Sigma_sims[t, k, , ]}, t)
    
    parms <- list()
    parms$beta_0 <- beta_0_temp
    parms$beta <- beta_temp
    parms$theta <- theta_temp
    parms$Sigma <- Sigma_temp
    parms$tau <- tau_temp
    parms$gamma <- gamma_temp
    
    z_star[i,] <- VARHMM_viterbi(obs, parms, pseudo = pseudo, L1_ball = L1_ball, 
                              draw = TRUE) 
    
    y_hat[i, , ] <- t(sapply(1:N, function(n) {rmvnorm(1, mu_sims[t, z_star[i,n], n, ], 
                                                       Sigma_temp[[z_star[i,n]]])}))
    
    log_pred_density[i,] <- t(sapply(1:N, function(n) {dmvnorm(obs[n,], mu_sims[t, z_star[i,n], n, ], 
                                                               Sigma_temp[[z_star[i,n]]], log = TRUE)}))
    
    i <- i + 1
  }
  
  z_hat <- VARHMM_viterbi(obs, parms_hat, pseudo = pseudo, L1_ball = L1_ball, 
                           draw = FALSE) 
  
  return(list(y_hat = y_hat, 
              z_hat = z_hat,
              z_star = z_star,
              log_pred_density = apply(log_pred_density, 2, logSumExp) - log(ndraw)))
}

VARHMM_get_predictive_sims <- function(sims, obs, pseudo = FALSE, L1_ball = TRUE, 
                                  ndraw = 50)
{
  #sims <- rstan::extract(fit)
  niter <- nrow(sims$beta)
  if(ndraw > niter) {stop("ndraw is larger than the number of  mcmc iterations:
                           choose a value between 1 and niter")}
  D <- ncol(obs)
  N <- nrow(obs)
  y_hat <- array(NA, dim = c(ndraw, N, D))
  log_pred_density <- array(NA, dim = c(ndraw, N))
  
  gamma_sims <- sims$gamma
  beta_0_sims <- sims$beta_0
  beta_sims <- sims$beta
  theta_sims <- get_theta_sims(sims, L1_ball)
  tau_sims <- sims$tau
  mu_sims <- get_mu_sims(sims, obs, L1_ball)
  Sigma_sims <- get_Sigma_sims(sims, pseudo)
  K <- dim(gamma_sims)[2]
  
  
  beta_0_hat <- lapply(1:K, function(kk) { 
    apply(beta_0_sims[, kk, ], 2, mean)})
  beta_hat <- lapply(1:K, function(kk) { 
    apply(beta_sims[, kk, ], 2, mean)})
  theta_hat <- lapply(1:K, function(kk) { 
    apply(theta_sims[, kk, ], 2, mean)})
  tau_hat <- lapply(1:K, function(kk) { 
    apply(tau_sims[, kk, ], 2, mean)})
  gamma_hat <- apply(gamma_sims,  c(2, 3), mean)
  Sigma_hat <- lapply(1:K, function(kk) {
    apply(Sigma_sims[, kk, , ], c(2, 3), mean)})
  
  
  parms_hat <- list()
  parms_hat$beta_0 <- beta_0_hat
  parms_hat$beta <- beta_hat
  parms_hat$theta <- theta_hat
  parms_hat$Sigma <- Sigma_hat
  parms_hat$tau <- tau_hat
  parms_hat$gamma <- gamma_hat
  
  
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
    gamma_temp <- gamma_sims[t, , ]
    Sigma_temp <- lapply(1:K, function(k, t) {Sigma_sims[t, k, , ]}, t)
    
    parms <- list()
    parms$beta_0 <- beta_0_temp
    parms$beta <- beta_temp
    parms$theta <- theta_temp
    parms$Sigma <- Sigma_temp
    parms$tau <- tau_temp
    parms$gamma <- gamma_temp
    
    z_star <- VARHMM_viterbi(obs, parms, pseudo = pseudo, L1_ball = L1_ball, 
                             draw = TRUE) 
    
    y_hat[i, , ] <- t(sapply(1:N, function(n) {rmvnorm(1, mu_sims[t, z_star[n], n, ], 
                                                       Sigma_temp[[z_star[n]]])}))
    
    log_pred_density[i,] <- t(sapply(1:N, function(n) {dmvnorm(obs[n,], mu_sims[t, z_star[n], n, ], 
                                                               Sigma_temp[[z_star[n]]], log = TRUE)}))
    
    i <- i + 1
  }
  
  z_hat <- VARHMM_viterbi(obs, parms_hat, pseudo = pseudo, L1_ball = L1_ball, 
                          draw = FALSE) 
  
  return(list(y_hat = y_hat, 
              z_hat = z_hat,
              log_pred_density = apply(log_pred_density, 2, logSumExp) - log(ndraw)))
}


# - 
quad_form_diag <- function(Omega, tau) {
  diag(tau) %*% Omega %*% diag(tau)
}


# - 
VARHMM_getBayesEstimates <- function(sims, pseudo = FALSE, L1_ball = TRUE)
{
  gamma_sims <- sims$gamma
  beta_0_sims <- sims$beta_0
  beta_sims <- sims$beta
  tau_sims <- sims$tau
  Sigma_sims <- get_Sigma_sims(sims, pseudo)
  theta_sims <- get_theta_sims(sims, L1_ball)
  K <- dim(gamma_sims)[2]
  
  
  beta_0_hat <- lapply(1:K, function(kk) { 
    apply(beta_0_sims[, kk, ], 2, mean)})
  beta_hat <- lapply(1:K, function(kk) { 
    apply(beta_sims[, kk, ], 2, mean)})
  tau_hat <- lapply(1:K, function(kk) { 
    apply(tau_sims[, kk, ], 2, mean)})
  gamma_hat <- apply(gamma_sims,  c(2, 3), mean)
  Sigma_hat <- lapply(1:K, function(kk) {
    apply(Sigma_sims[, kk, , ], c(2, 3), mean)})
  
  
  parms_hat <- list()
  parms_hat$beta_0 <- beta_0_hat
  parms_hat$beta <- beta_hat
  parms_hat$Sigma <- Sigma_hat
  parms_hat$tau <- tau_hat
  parms_hat$gamma <- gamma_hat
  
  
  if (L1_ball) {
    theta_hat <- lapply(1:K, function(kk) { 
      apply(theta_sims[, kk, ], 2, mean)})
    parms_hat$theta <- theta_hat
  }
  
  return(parms_hat)
}

# - 
VARHMM_forwardBackwards <- function(obs, allprobs, gamma)
{
  K <- nrow(gamma)
  N <- nrow(obs)
  
  lalpha <- lbeta <- matrix(NA, K, N)
  gamma_init <- solve(t(diag(K)-gamma+1),rep(1,K))
  foo <- gamma_init*allprobs[1, ]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[, 1] <- log(foo) + lscale
  
  for (i in 2:N) {
    foo <- foo%*%gamma*allprobs[i,]
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
    lalpha[, i] <- log(foo) +lscale
  }
  llk <- lscale
  
  lbeta[, N] <- rep(0, K)
  foo <- rep(1/K, K)
  lscale <- log(K)
  
  for (i in (N-1):1) {
    foo <- gamma%*%(allprobs[i+1,]*foo)
    lbeta[,i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  return(list(la = lalpha, lb = lbeta, llk = llk))
}

# - 
VARHMM_stateProbs <- function(sims, obs, pseudo = FALSE, L1_ball = TRUE)
{
  parms_hat <- VARHMM_getBayesEstimates(sims, pseudo, L1_ball)
  K <- nrow(parms_hat$gamma)
  N <- nrow(obs)
  
  allprobs <- VAR_emissions(obs, parms_hat, pseudo, L1_ball)
  fb <- VARHMM_forwardBackwards(obs, allprobs, parms_hat$gamma)
  la <- fb$la
  lb <- fb$lb
  llk <- fb$llk
  
  stateprobs <- matrix(NA, ncol = N, nrow = K)
  for (n in 1:N) stateprobs[, n] <- exp(la[, n] + lb[, n] - llk)
  return(stateprobs)
}

# The previous version of this function just does so for the posterior means, 
# why not use the full posterior sample 
# - 
VARHMM_stateProbs_Bayesian <- function(sims, obs, pseudo = FALSE, L1_ball = TRUE, ndraw = 50)
{
  #sims <- rstan::extract(fit)
  niter <- nrow(sims$beta)
  if(ndraw > niter) {stop("ndraw is larger than the number of  mcmc iterations:
                           choose a value between 1 and niter")}
  gamma_sims <- sims$gamma
  beta_0_sims <- sims$beta_0
  beta_sims <- sims$beta
  tau_sims <- sims$tau
  Sigma_sims <- get_Sigma_sims(sims, pseudo)
  theta_sims <- get_theta_sims(sims, L1_ball)
  
  N_MCMC <- length(sims$lp__)
  K <- dim(gamma_sims)[2]
  N <- nrow(obs)
  stateprobs <- array(NA, dim = c(ndraw, K, N))
  
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
    parms_hat$gamma <- gamma_sims[t,,]
    if (L1_ball) {
      parms_hat$theta <- lapply(1:K, function(kk) {theta_sims[t,kk,]}) 
    }
      
    allprobs <- VAR_emissions(obs, parms_hat, pseudo, L1_ball)
    fb <- VARHMM_forwardBackwards(obs, allprobs, parms_hat$gamma)
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
VARHMM_plotStateProbs <- function(fit, obs, pseudo = FALSE, L1_ball = TRUE)
{
  sims <- rstan::extract(fit)
  K <- dim(sims$gamma)[2]
  N <- nrow(obs)
  time <- 1:N
  if (K>2) {
    z_col <- brewer.pal(n = K, name = "Spectral")
  } else {
    z_col <- c("lightsalmon", "lightblue1")
  }
  
  state_probs <- VARHMM_stateProbs(sims, obs, pseudo, L1_ball)
  parms_hat <- VARHMM_getBayesEstimates(sims, pseudo, L1_ball)
  z_hat <- VARHMM_viterbi(obs, parms_hat, pseudo, L1_ball, draw = FALSE) 
  z_aux <- sapply(1:N, function(n) z_col[z_hat[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
  
  plot(1:N, type = "n", xlab = "Time", ylab = "Prob state", 
       cex = 0.5, las = 1,
       xlim = c(1, N),
       ylim = c(0,1), cex.lab = 1.1)
  
  
  plot_p <- matrix(NA,nrow=N-1,ncol=2*K)
  a <- 0
  for (n in 2:N) {
    
    for (j in 1:K) {
      plot_p[n-1, (j*2-1)] <- state_probs[j, n-1]
      plot_p[n-1, j*2] <- state_probs[j, n]
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



