c_light = "gray95"
c_light_highlight = "gray90"
c_mid = "gray85"
c_mid_highlight = "gray80"
c_dark = "gray75"


# ? 
get_companion_matrix <- function(A)
{
  P <- dim(A)[3]
  D <- dim(A)[1]
  
  A_aux <- rep(0, D)
  for (p in 1:P) {
    A_aux <- cbind(A_aux, A[, , p])
  }
  A_aux <- A_aux[, -1]
  I_aux <- cbind(do.call("cbind", replicate(P-1, diag(D), simplify=FALSE)), 
                 matrix(0, D, D))
  A_temp <- do.call("rbind", replicate(P-1, I_aux, simplify=FALSE))
  A_companion <- rbind(A_aux, A_temp)
  
  return(A_companion)
}

#- ? 
generate_sparse_beta <- function(D, P, sparsity_level, 
                                 abs_beta_low = 0.2, abs_beta_up = 0.8) 
{
  beta <- array(NA, dim = c(D, D, P))
  #zeros <- sample(c(0, 1), P*D^2, replace = TRUE,
  #                prob = c(1 - sparsity_level, sparsity_level))
  stab <- FALSE
  
  while (stab == FALSE) {
    for (p in 1:P) {
      zeros <- sample(c(0, 1), D^2, replace = TRUE,
                      prob = c(1 - sparsity_level, sparsity_level))
      sign <- (2*rbinom(D^2, size = 1, prob = 0.5) - 1)
      elements <- (1 - zeros) * sign * runif(D^2, abs_beta_low, abs_beta_up)
      beta[, , p] <-  matrix(elements, D, D, byrow = T)
    }
    A_comp <- get_companion_matrix(beta)
    stab <- all(abs(eigen(A_comp)$values) < 1)
  }
  
  return(beta)
}

# - 
generate_covariance <- function(D) {
  p <- qr.Q(qr(matrix(rnorm(D^2), D)))
  Sigma <- crossprod(p, p * sample(1:3, D, replace = T))
  return(Sigma)
}



# - 
quad_form_diag <- function(Omega, tau) {
  diag(tau) %*% Omega %*% diag(tau)
}


# - 
get_Sigma_sims <- function(sims, pseudo = FALSE) 
{
  tau_sims <- sims$tau
  n_MCMC <- dim(tau_sims)[1]
  K <-  dim(tau_sims)[2]
  D <- dim(tau_sims)[3]
  
  Sigma_sims <- array(NA, dim = c(n_MCMC, K, D, D))
  
  if (pseudo) {
    for (tt in 1:n_MCMC) {
      tau_all <- tau_sims[tt, , ]
      for (kk in 1:K) {
        tau <- tau_all[kk, ]
        Sigma_sims[tt, kk, , ] <- diag(tau) # tau is modeled as variances
      }
    }
  } else {
    for (kk in 1:K) {
      for (tt in 1:n_MCMC) {
        Omega <- sims$L_Omega[tt, kk, ,] %*% t(sims$L_Omega[tt, kk, ,])
        tau <- sims$tau[tt, kk, ] # tau is modeled as scale
        Sigma_sims[tt, kk, , ] <- quad_form_diag(Omega, tau)
      }
    }
  }
  return(Sigma_sims)
}

# 
get_mu_sims <- function(sims, obs, L1_ball = FALSE)
{
  
  beta_0_sims <- sims$beta_0  
  beta_sims <- sims$beta
  tau_sims <- sims$tau
  theta_sims <- sims$theta  
  
  K <- dim(sims$beta_0)[2]
  D <- ncol(obs)
  N <- nrow(obs)
  P_hat <- dim(beta_sims)[3]
  P <- P_hat/(D*D)
  n_MCMC <- dim(beta_sims)[1]
  
  mu_sims <- array(NA, dim = c(n_MCMC, K, N, D))
  
  for (tt in 1:n_MCMC) {
    for (kk in 1:K) {
      for (n in 1:P) {
        mu_sims[tt, kk, n, ] <- obs[n, ]
      }
      
      for(n in (P+1):N) {
        mu <- beta_0_sims[tt, kk, ]
        for(p in 1:P) {
          if (L1_ball) {
            beta_temp <- get_beta_temp(theta_sims[tt, kk,  ], p, D)
          } else {
            beta_temp <- get_beta_temp(beta_sims[tt, kk,  ], p, D)
          }
          mu <- mu + beta_temp %*% obs[n-p, ]
        }
        mu_sims[tt, kk, n, ] <- mu
      }
    }
  }
  
  return(mu_sims)
}

# - input beta is a D*D*P vector
get_beta_temp <- function(beta, p, D)
{
  a <- ifelse(p == 1, 1, ((D*D)*(p-1)) + 1)
  b <- (D*D)*p
  return(matrix(beta[a:b], D, D, byrow = T))
}

# - 
standardize <- function(obs)
{
  N <- nrow(obs)
  D <- ncol(obs)
  out <- matrix(NA, nrow = N, ncol = D)
  for (dd in 1:D) {
    out[, dd] <- (obs[, dd] - mean(obs[, dd]))/sd(obs[, dd])
  }
  return(out)
}


# - 
l1_ball_projection <- function(beta, r){
  p <- length(beta)
  abs_beta <- abs(beta)
  abs_beta_sort_index <- order(abs_beta, decreasing = TRUE)
  mu <- pmax(cumsum(abs_beta[abs_beta_sort_index]) - r, 0)
  c <- max(which((abs_beta[abs_beta_sort_index] > mu/(1:p)) == TRUE))
  theta <- sign(beta)*pmax(abs(beta) - mu[c]/c, 0)
  return(theta)
}

# - 
get_theta_sims <- function(sims, L1_ball = FALSE, tilde_mu = FALSE)##CHANGE THIS BACK
{
  
  N_MCMC <- dim(sims$beta_0)[1]
  K <- dim(sims$beta_0)[2]
  D <- dim(sims$beta_0)[3]
  P <- (dim(sims$beta)[3])/(D*D)
  beta_sims <- sims$beta
  
  if (L1_ball) {
    if (!tilde_mu){
      ray_sims <- sims$ray
      theta_sims <- array(NA, dim = c(N_MCMC, K, D*D*P))
      for (kk in 1:K) {
        for (tt in 1:N_MCMC) {
          theta_sims[tt, kk, ] = l1_ball_projection(beta_sims[tt, kk, ], ray_sims[tt, kk])
        }
      }
    }
    if (tilde_mu){
      tilde_mu_sims <- sims$tilde_mu
      theta_sims <- array(NA, dim = c(N_MCMC, K, D*D*P))
      for (kk in 1:K) {
        for (tt in 1:N_MCMC) {
          theta_sims[tt, kk, ] = l1_ball_projection2(beta_sims[tt, kk, ], tilde_mu_sims[tt, kk])
        }
      }
    }
    
    
    return(theta_sims)
  } else {
    return(beta_sims)
  }
}

# - 
get_beta_est <- function(sims, L1_ball = FALSE, mat = TRUE)
{
  K <- dim(sims$beta_0)[2]
  D <- dim(sims$beta_0)[3]
  P <- (dim(sims$beta)[3])/(D*D)
  
  theta_sims <- get_theta_sims(sims, L1_ball)
  theta_est_all <- apply(theta_sims, c(2, 3), mean)
  
  theta_mat <- array(NA, c(K, D, D, P))
  
  for (kk in 1:K) {
    for (p in 1:P) {
      theta_mat[kk, , , p] = get_beta_temp(theta_est_all[kk, ], p, D)
    }
  }
  
  if (mat) {
    return(theta_mat)
  } else {
    return(theta_est_all)
  }
}


# - ? 
get_inclusion_sims <- function(sims, p_est = TRUE)
{
  theta_sims_temp <- get_theta_sims(sims, L1_ball = TRUE)
  
  N_MCMC <- dim(sims$beta_0)[1]
  K <- dim(sims$beta_0)[2]
  D <- dim(sims$beta_0)[3]
  P <- (dim(sims$beta)[3])/(D*D)
  
  inclusion <- array(NA, dim = c(N_MCMC, K, D, D, P))
  
  for (tt in 1:N_MCMC) {
    for (kk in 1:K) {
      for (pp in 1:P) {
        inclusion[tt, kk, , , pp] <- get_beta_temp(theta_sims_temp[tt, kk, ] == 0,
                                                   pp, D)
      }
    }
  }
  
  if (p_est) {
    p_est <- array(NA, dim = c(K, D, D, P))
    for (kk in 1:K) {
      for (pp in 1:P) {
        p_est[kk, , , pp] = 1 - apply(inclusion[, kk, ,  , pp], c(2, 3), mean)
      }
    }
  }
  
  
  return(list(inclusion = inclusion, p_est = p_est))
}


# 
get_mu_sims <- function(sims, obs, L1_ball = FALSE)
{
  
  beta_0_sims <- sims$beta_0  
  beta_sims <- sims$beta
  tau_sims <- sims$tau
  theta_sims <- get_theta_sims(sims, L1_ball) 
  

  K <- dim(sims$beta_0)[2]
  D <- ncol(obs)
  N <- nrow(obs)
  P_hat <- dim(beta_sims)[3]
  P <- P_hat/(D*D)
  n_MCMC <- dim(beta_sims)[1]
  
  mu_sims <- array(NA, dim = c(n_MCMC, K, N, D))
  
  for (tt in 1:n_MCMC) {
    for (kk in 1:K) {
      for (n in 1:P) {
        mu_sims[tt, kk, n, ] <- obs[n, ]
      }
      
      for(n in (P+1):N) {
        mu <- beta_0_sims[tt, kk, ]
        for(p in 1:P) {
          if (L1_ball) {
            beta_temp <- get_beta_temp(theta_sims[tt, kk,  ], p, D)
          } else {
            beta_temp <- get_beta_temp(beta_sims[tt, kk,  ], p, D)
          }
          mu <- mu + beta_temp %*% obs[n-p, ]
        }
        mu_sims[tt, kk, n, ] <- mu
      }
    }
  }
  
  return(mu_sims)
}

#  simplex to tpm 
a_to_mat <- function(a) 
{
  if (length(a) == 2) {
    K <- 2
  } else {
    K <- nrow(a)
  }
  gamma_out <- matrix(0, K, K)
  gamma_out[!diag(K)] <- as.vector(t(a))
  t(gamma_out)
  return(t(gamma_out))
}

# get gamma sims from a_sims
get_gamma_sims <- function(a) 
{
  n_iter = nrow(a)
  get.gamma <- function(t, a) {a_to_mat(a[t, , ])}
  simplify2array(lapply(1:n_iter, get.gamma, a = a))
}


# - ?
VAR_emissions <- function(obs, parms, pseudo = FALSE, L1_ball = TRUE)
{
  beta_0 <- parms$beta_0
  beta <- parms$beta
  theta <- parms$theta
  tau <- parms$tau
  K <- length(beta_0)
  
  
  if (pseudo) {
    Sigma <- lapply(1:K, function(kk) {diag(tau[[kk]])})
  } else {
    Sigma <- parms$Sigma
  }
  
  N <- nrow(obs)
  D <- ncol(obs)
  P_hat <- length(beta[[1]])
  P <- P_hat/(D*D)
  
  allprobs <- matrix(NA, N, K)
  
  # emission first P obs
  for (k in 1:K) {
    for (n in 1:P) {
      allprobs[n, k] = dmvnorm(obs[n, ], obs[n, ], Sigma[[k]])
    }
  }
  
  # emissions (P+1):N
  for (k in 1:K) {
    for (n in (P+1):N) {
      mu = beta_0[[k]]
      for (p in 1:P) {
        if (L1_ball) {
          beta_temp = get_beta_temp(theta[[k]], p, D);
        } else {
          beta_temp = get_beta_temp(beta[[k]], p, D);
        }
        mu = mu + beta_temp %*% obs[n - p, ];
      }
      allprobs[n, k] = dmvnorm(obs[n, ], mu, Sigma[[k]])
    }
  }
  
  return(allprobs)
  
}


# - 
plotPosteriorPredictive <- function(obs, y_hat, z_hat, K)
{
  
  
  ndraw = dim(y_hat)[1]
  D <- dim(y_hat)[3]
  N <- dim(y_hat)[2]
  
  
  if (D == 8) {
    n_row = 4
    n_col = 2
  } else {
    n_row = 5
    n_col = ceiling(D/n_row)
    par(mfrow = c(n_row, n_col))
  }
  
  # 
  # n_row = 5
  # n_col = ceiling(D/n_row)
  # par(mfrow = c(n_row, n_col))
  
  if (K>2) {
    z_col <- brewer.pal(n = K, name = "Spectral")
  } else {
    z_col <- c("lightsalmon", "lightblue1")
  }
  
  z_aux <- sapply(1:N, function(n) z_col[z_hat[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
  col_all <- rainbow(D)
  
  par(mfrow = c(n_row, n_col), 
      mai = c(0.4, 0.55, 0.1, 0.3))
  
  for (dd in 1:D) {
    
    obs_univ <- obs[, dd]
    plot.new()
    plot.window(ylim = c(min(obs_univ) - 0.6, max(obs_univ) + 0.6), 
                xlim = c(1, N),  xlab = "Time",  ylab = paste("y_", dd, sep=""), 
                cex.lab = 1.1)
    title(xlab = "Time", ylab = paste("y_", dd, sep=""))
    for (j in 1:(length(cp_loc) - 1)) { 
      rect(cp_loc[j], ceiling(min(obs_univ) - 2), cp_loc[j+1], ceiling(max(obs_univ) + 1), 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
    }
    
    probs <- seq(from=0.1, to=0.9, by=0.1)
    
    t = 1
    
    cred <- sapply(1:N, function(t) quantile(y_hat[, t, dd], probs = probs))
    polygon(c(1:N, rev(1:N)), c(cred[1,], rev(cred[9,])),
            col = c_light, border = NA)
    polygon(c(1:N, rev(1:N)), c(cred[2,], rev(cred[8,])),
            col = c_light_highlight, border = NA)
    polygon(c(1:N, rev(1:N)), c(cred[3,], rev(cred[7,])),
            col = c_mid, border = NA)
    polygon(c(1:N, rev(1:N)), c(cred[4,], rev(cred[6,])),
            col = c_mid_highlight, border = NA)
    lines(1:N, cred[5,], col=c_dark, lwd=1)
    points(1:N, obs_univ, pch = 16, cex = 0.9, col = "black")
    points(1:N, obs_univ, pch = 16, cex = 0.8, col = z_aux)
    axis(1); axis(2)
  }
}

# - 
plotDAG <- function(X, ylabels, color, main) {
  
  X[X<0.5] = 0
  X[X>=0.5] = 1
  D <- nrow(X)
  
  X_star <- matrix(0, nrow = D*2, ncol = D*2)
  for(ii in 1:D) {
    for (jj in 1:D) {
      X_star[2*ii-1, 2*jj] = t(X)[ii, jj] 
    }
  }
  graph_mat<- as(X_star, "sparseMatrix"); graph_mat
  label_graph <-  vector(mode='character',length=0)
  for (dd in 1:D) {
    label_graph <- c(label_graph, paste(ylabels[dd], "(t-1)", sep = ""))
    label_graph <-  c(label_graph, paste(ylabels[dd], "(t)", sep = ""))
  }
  dimnames(graph_mat) = list(label_graph, label_graph)
  
  G <- graph_from_adjacency_matrix(
    graph_mat,
    mode = "directed"
  )
  G_l <- layout_on_grid(G, height = D, width = 2, dim = 2)
  plot(G,  layout = G_l, vertex.size = 36.5, 
       edge.arrow.size = 1.3, 
       vertex.label.cex = 1.1, 
       margin = 0.1, 
       vertex.color = color)
  # "lightsalmon"
  # "lightblue1"
  title(main, cex.main = 4)
}



# - 
plotPosteriorCorrelation <- function(corr_sims, z_hat, K, corr_label)
{
  
  N <- length(z_hat)
  if (K>2) {
    z_col <- brewer.pal(n = K, name = "Spectral")
  } else {
    z_col <- c("lightblue1", "lightsalmon")
  }
  
  z_aux <- sapply(1:N, function(n) z_col[z_hat[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
  
  ylabels <- paste("Corr", "(", corr_label[1], " / ", 
                   corr_label[2], ")", sep = "")
  
  # careful this is specific to our application, correlation can be negative!!
  my_ylimit <- c(-0.1, 1)
  
  plot(1:N, type = "n",  ylab = ylabels, 
       xlim = c(1, N),
       ylim = my_ylimit, cex.lab = 1.1, xlab = "")
  
  probs <- seq(from=0.1, to=0.9, by=0.1)
  cred <- sapply(1:N, function(t) quantile(corr_sims[, t], probs = probs))
  
  # for (j in 1:(length(cp_loc) - 1)) { 
  #   rect(cp_loc[j], ceiling(min(corr_sims) - 2), cp_loc[j+1], ceiling(max(corr_sims) + 1), 
  #        col =scales::alpha(z_aux[cp_loc[j+1]], 0.2), border = F)
  # }
  for (j in 1:(length(cp_loc) - 1)) { 
    rect(cp_loc[j], my_ylimit[1], cp_loc[j+1], my_ylimit[2], 
         col =scales::alpha(z_aux[cp_loc[j+1]], 0.2), border = F)
  }
  
  
  
  polygon(c(1:N, rev(1:N)), c(cred[1,], rev(cred[9,])),
          col = scales::alpha(c_light, 0.8), border = NA)
  polygon(c(1:N, rev(1:N)), c(cred[2,], rev(cred[8,])),
          col = scales::alpha(c_light_highlight, 0.8), border = NA)
  polygon(c(1:N, rev(1:N)), c(cred[3,], rev(cred[7,])),
          col = scales::alpha(c_mid, 0.8), border = NA)
  polygon(c(1:N, rev(1:N)), c(cred[4,], rev(cred[6,])),
          col = scales::alpha(c_mid_highlight, 0.9), border = NA)
  lines(1:N, cred[5,], col= scales::alpha(c_dark, 1.0), lwd=3)
}



#' Alejandra Avalos Pacheco's function to splot sparsity 

#' Generate heatmap of matrices
#' @param Matrix matrix to plot
#' @param limit limit for the values of the plot
#' @export
#' @example
#' library(ggplot2)
#' sigma = matrix(rbinom(10,1,.30), 5, 2)
#' plot.heat(sigma,limit=c(-1,1))
plot.heat <- function(Matrix,Xlab="",Ylab="",Main="",limit=c(-2,2)){
  Matrix = as.matrix(Matrix)
  colnames(Matrix)<-NULL
  rownames(Matrix)<-NULL
  x = reshape2::melt(data.frame(Matrix,ind=c(nrow(Matrix):1)),id="ind")
  colnames(x)<-c("X1","X2","value")
  p_X_heat = ggplot(data = x, aes(x=X2, y=X1, fill=value)) +
    theme_bw() +
    #geom_tile(show.legend = F) +
    geom_tile() +
    xlab(Xlab) +
    ylab(Ylab) +
    ggtitle(Main) +
    theme(axis.title=element_text(size=14,face="bold"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    scale_fill_gradient2(limits=limit) + 
    theme(legend.position="bottom",
          legend.margin=margin(t=-25))
  return(p_X_heat)
}

# Unstandardising the beta coefficients of a VAR fitted to standardised data
beta_unstandardise <- function(beta, bar_mu, bar_s){
  # beta is D x D x P
  # bar_mu is a D dimensional vector of dimension means
  # bar_s is a D dimensional vector of dimension s.d.
  dim_beta <- dim(beta)
  D <- dim_beta[1]
  P <- dim_beta[3]
  if(is.na(P) == TRUE){
    P <- 1
    beta <- array(beta, dim = c(D, D, P))
  }
  
  beta_unsd <- array(NA, dim = c(D, D, P))
  for(d in 1:D){
    for(j in 1:D){
      for(p in 1:P){
        beta_unsd[d, j, p] <- bar_s[d]*beta[d, j, p]/bar_s[j]
      }
    }
  }
  return(beta_unsd)
}

simulation_plot <- function(data, plt_all = TRUE, plt_univ = FALSE)
{
  obs <- data$obs
  D <- ncol(obs)
  
  state_seq <- data$state
  state_phase <- data$state_phase
  K <- length(unique(state_seq))
  
  K <- length(unique(state_seq))
  if (K>2) {
    z_col <- brewer.pal(n = K, name = "Spectral")
  } else {
    z_col <- c( "lightblue1", "lightsalmon")
  }
  
  z_aux <- sapply(1:N, function(n) z_col[state_seq[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
  col_all <- rainbow(D)
  
  if (plt_all) {
    
    
    # - multivariate single plot
    par(mfrow = c(1, 1))
    title_main <- ""
    plot(1:N, type = "n",  ylab =  "y", xlab = "Time",
         cex = 0.5, las = 1,
         xlim = c(1, N),
         ylim = c(min(obs) - 1, max(obs) + 1), cex.lab = 1.1, 
         main = title_main)
    
    for (j in 1:(length(cp_loc) - 1)) { 
      rect(cp_loc[j], ceiling(min(obs) - 2) , cp_loc[j+1], ceiling(max(obs) + 2) , 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.4), border = F)
    }
    for (d in 1:D) {
      lines(obs[, d], pch = 20,
            col = col_all[d], cex = 0.6, type = "o")
    }
    
    legend_aux <- numeric(K)
    for (kk in 1:K) {
      legend_aux[kk] <- paste(unique(state_seq)[kk], 
                              unique(state_phase)[kk], sep = ": ")
    }
    legend("topright", 
           legend  = legend_aux, 
           pch = 15, 
           col = unique(z_aux), 
           cex = 0.5,
           bg = "white")
  }
  
  if (plt_univ) {
    #readline(prompt="Plot, univariate: Press [enter] to continue")
    
    if (only_acc) {
      n_row = 4
      n_col = 2
    } else {
      n_row = 5
      n_col = ceiling(D/n_row)
    }
    par(mfrow = c(n_row, n_col))
    
    par(mfrow = c(n_row, n_col), 
        mai = c(0.4, 0.55, 0.1, 0.3))
    
    for (dd in 1:D) {
      obs_univ <- obs[, dd]
      plot.new()
      plot.window(ylim = c(min(obs_univ) - 0.3, max(obs_univ) + 0.3), 
                  xlim = c(1, N),  xlab = "Time",  ylab = paste("y_", dd, sep=""), 
                  cex.lab = 1.1)
      #y_name = names(data_gesture)[1:D][dd]
      # y_name = rep("", D) # need to fix this
      title(xlab = "Time", ylab = paste("y_", dd,sep=""))
      for (j in 1:(length(cp_loc) - 1)) { 
        rect(cp_loc[j], ceiling(min(obs_univ)-2), cp_loc[j+1], ceiling(max(obs_univ) + 2), 
             col =scales::alpha(z_aux[cp_loc[j+1]], 0.2), border = F)
      }
      lines(1:N, obs_univ, pch = 20, cex = 0.6, col = col_all[dd], type = "o")
      axis(1); axis(2)
      
      legend_aux <- numeric(K)
      for (kk in 1:K) {
        legend_aux[kk] <- paste(unique(state_seq)[kk], 
                                unique(state_phase)[kk], sep = ": ")
      }
      legend("topright", 
             legend  = legend_aux, 
             pch = 15, 
             col = unique(z_aux), 
             cex = 0.5)
    }
  }
  
  
  
  #return()
}