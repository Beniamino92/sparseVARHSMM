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

# causal filter
filter_causal <- function(y, win_s) {
  
  N <- length(y)
  N_smooth <- N - win_s + 1
  y_smooth <- c()
  for (t in 1:N_smooth) {
    y_smooth <- c(y_smooth, sum(y[t:(t+win_s-1)])/win_s)
  }
  y_smooth <- c(rep(NA, win_s  - 1), y_smooth)
  
  return(y_smooth)
}

# - get data Gesture Phase

get_data_Gesture_Phase <- function(path, patient, story, processed, 
                                   down_sample = 1, win_s = 0, 
                                   rest_action = TRUE,
                                   scale = TRUE, sqrt_transf = FALSE,
                                   plt_all = FALSE, plt_univ = FALSE,
                                   plt_only_acc = TRUE)
{
  path_story <- paste(patient, story, 
                      ifelse(processed, "_va3.csv", "_raw.csv"), 
                      sep = "")
  path_data <- paste(path, path_story, sep = "")
  
  data_gesture <- read.csv(path_data)
  N_all <- dim(data_gesture)[1]
  D <- ifelse(processed, dim(data_gesture)[2] - 1, dim(data_gesture)[2] - 2)
  
  obs <- data_gesture[, 1:D]
  state_seq <- as.integer(factor(data_gesture$phase))
  state_phase <- data_gesture$phase
  N <- length(state_seq)
  
  if (processed) {
    state_seq <- as.integer(factor(data_gesture$Phase))
    state_phase <- data_gesture$Phase
    N <- length(state_seq)
    if (sqrt_transf) {
      for (dd in (D-7):D) {
        obs[, dd] <- sqrt(obs[, dd])
      }
    }
  }
  
  # only rest and action
  if (rest_action) {
    state_seq[which(state_phase == "D")] = 1
    state_seq[which(state_phase != "D")] = 2
      
    state_phase[which(state_phase == "D")] = "Rest"
    state_phase[which(state_phase != "Rest")] = "Active"
  }
  
  
  # smoothing
  if (win_s > 0) {
    N <- length(state_seq)
    out <- matrix(NA, nrow = N, ncol = D)
    for (dd in 1:D) {
      out[, dd] <- filter_causal(obs[, dd], win_s)
    }
    obs <- out[win_s:N, ]
    state_seq = state_seq[win_s:N]
    N <- length(state_seq)
  }
  
  # standardize
  if (scale) {
    obs <- standardize(obs)
  } 
  
  # sub-sampling 
  idxs <- seq(1, N, by = down_sample)
  obs <- obs[idxs, ]
  state_seq <- state_seq[idxs]
  state_phase <- state_phase[idxs]
  N <- length(state_seq)
  
  
  # preparing for plot 
  K <- length(unique(state_seq))
  if (K>2) {
    z_col <- brewer.pal(n = K, name = "Accent")
  } else {
    z_col <-c( "lightblue1", "lightsalmon")
  }
  col_all <- rainbow(D)
  z_aux <- sapply(1:N, function(n) z_col[state_seq[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
  col_all <- rainbow(D)

  
  
  
  if (plt_only_acc) {
    
    col_all <- rainbow(n = 8)
    n_row = 4
    n_col = 2
    par(mfrow = c(n_row, n_col))
    
    par(mfrow = c(n_row, n_col), 
        mai = c(0.4, 0.55, 0.1, 0.3))
    
    for (dd in (D-7):D) {
      obs_univ <- obs[, dd]
      plot.new()
      plot.window(ylim = c(min(obs_univ) - 0.3, max(obs_univ) + 0.3), 
                  xlim = c(1, N),  xlab = "Time",  ylab = paste("y_", dd-24, sep=""), 
                  cex.lab = 1.1)
      y_name = names(data_gesture)[1:D][dd]
      title(xlab = "Time", ylab = paste("y_", dd-24, sep=""))
      for (j in 1:(length(cp_loc) - 1)) { 
        rect(cp_loc[j], ceiling(min(obs_univ)-2), cp_loc[j+1], ceiling(max(obs_univ) + 2), 
             col =scales::alpha(z_aux[cp_loc[j+1]], 0.2), border = F)
      }
      lines(1:N, obs_univ, pch = 20, cex = 0.6, col = col_all[dd-24], type = "o")
      axis(1); axis(2)
      legend_aux <- numeric(K)
      for (kk in 1:K) {
        legend_aux[kk] <- paste(unique(state_seq)[kk], 
                                unique(state_phase[idxs])[kk], sep = ": ")
      }
      for (j in 1:(length(cp_loc) - 1)) { 
        rect(cp_loc[j], ceiling(min(obs) - 2) , cp_loc[j+1], ceiling(max(obs) + 2) , 
             col =scales::alpha(z_aux[cp_loc[j+1]], 0.4), border = F)
      }
      legend("topright", 
             legend  = legend_aux, 
             pch = 15, 
             col = unique(z_aux), 
             cex = 0.5)
    }
    
    
    
    
    readline(prompt="Plot, univariate: Press [enter] to continue")
    dev.off()
    par(mfrow = c(1, 1))
    title_main <- paste("Patient ", patient, "; Story ", story, sep = "")
    plot(1:N, type = "n",  ylab = "Sqrt(Physical Activity)", xlab = "Time",
         cex = 0.5, las = 1,
         xlim = c(1, N),
         ylim = c(min(obs[, (D-7):D]) - 1, max(obs[, (D-7):D]) + 1), cex.lab = 1.1, 
         main = title_main)
    
    for (d in (D-7):D) {
      lines(obs[, d], pch = 20,
            col = col_all[d-24], cex = 0.6, type = "o")
    }
    for (j in 1:(length(cp_loc) - 1)) { 
      rect(cp_loc[j], ceiling(min(obs) - 2) , cp_loc[j+1], ceiling(max(obs) + 2) , 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.4), border = F)
    }
    legend_aux <- numeric(K)
    for (kk in 1:K) {
      legend_aux[kk] <- paste(unique(state_seq)[kk], 
                              unique(state_phase[idxs])[kk], sep = ": ")
    }
    legend("topright", 
           legend  = legend_aux, 
           pch = 15, 
           col = unique(z_aux), 
           cex = 0.8)
    
  } else {
    
    if (plt_all) {
      
      # - multivariate single plot
      par(mfrow = c(1, 1))
      title_main <- paste("Patient ", patient, "; Story ", story, sep = "")
      plot(1:N, type = "n",  ylab = "", xlab = "Time",
           cex = 0.5, las = 1,
           xlim = c(1, N),
           ylim = c(min(obs) - 1, max(obs) + 1), cex.lab = 1.1, 
           main = title_main)
      
      for (d in 1:D) {
        lines(obs[, d], pch = 20,
              col = col_all[d], cex = 0.6, type = "o")
      }
      for (j in 1:(length(cp_loc) - 1)) { 
        rect(cp_loc[j], ceiling(min(obs) - 2) , cp_loc[j+1], ceiling(max(obs) + 2) , 
             col =scales::alpha(z_aux[cp_loc[j+1]], 0.4), border = F)
      }
      legend_aux <- numeric(K)
      for (kk in 1:K) {
        legend_aux[kk] <- paste(unique(state_seq)[kk], 
                                unique(state_phase[idxs])[kk], sep = ": ")
      }
      legend("topright", 
             legend  = legend_aux, 
             pch = 15, 
             col = unique(z_aux), 
             cex = 0.8)
    }
    
    if (plt_univ) {
      
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
        y_name = names(data_gesture)[1:D][dd]
        title(xlab = "Time", ylab = paste("y_", dd, ": ", y_name, sep=""))
        for (j in 1:(length(cp_loc) - 1)) { 
          rect(cp_loc[j], ceiling(min(obs_univ)-2), cp_loc[j+1], ceiling(max(obs_univ) + 2), 
               col =scales::alpha(z_aux[cp_loc[j+1]], 0.2), border = F)
        }
        lines(1:N, obs_univ, pch = 20, cex = 0.2, col = col_all[dd], type = "o")
        axis(1); axis(2)
      }
    }
  }
  
  return(list(obs = obs, 
              state = state_seq, 
              state_phase = state_phase))
}


# - subset data Gesture Phase, given inteval 
subset_data_Gesture_Phase <- function(Gesture_Phase, interval, only_acc = TRUE, 
                                      plt_all = TRUE, plt_univ = FALSE)
{
  obs <- Gesture_Phase$obs
  D <- ncol(obs)
  
  if (only_acc) {obs <- obs[, (D-7):D]}
  
  state_seq <- Gesture_Phase$state
  state_phase <- Gesture_Phase$state_phase
  K <- length(unique(state_seq))
  
  
  obs_sub <- obs[interval[1]:interval[2], ]
  state_seq_sub <- state_seq[interval[1]:interval[2]]
  state_phase_sub <- state_phase[interval[1]:interval[2]]
  N_sub <- nrow(obs_sub)
  D <- ncol(obs_sub)
  
  if (plt_all) {
    
    K <- length(unique(state_seq))
    if (K>2) {
      z_col <- brewer.pal(n = K, name = "Accent")
    } else {
      z_col <- c( "lightblue1", "lightsalmon")
    }
    # - multivariate single plot
    par(mfrow = c(1, 1))
    title_main <- paste("Patient ", patient, "; Story ", story, sep = "")
    plot(1:N_sub, type = "n",  ylab =  "Sqrt(Physical Activity)", xlab = "Time",
         cex = 0.5, las = 1,
         xlim = c(1, N_sub),
         ylim = c(min(obs_sub) - 1, max(obs_sub) + 1), cex.lab = 1.1, 
         main = title_main)
    
    z_aux <- sapply(1:N_sub, function(n) z_col[state_seq_sub[n]])
    cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N_sub)
    col_all <- rainbow(D)
    for (j in 1:(length(cp_loc) - 1)) { 
      rect(cp_loc[j], ceiling(min(obs_sub) - 2) , cp_loc[j+1], ceiling(max(obs_sub) + 2) , 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.4), border = F)
    }
    for (d in 1:D) {
      lines(obs_sub[, d], pch = 20,
            col = col_all[d], cex = 0.6, type = "o")
    }
    
    legend_aux <- numeric(K)
    for (kk in 1:K) {
      legend_aux[kk] <- paste(unique(state_seq_sub)[kk], 
                              unique(state_phase)[kk], sep = ": ")
    }
    legend("topright", 
           legend  = legend_aux, 
           pch = 15, 
           col = unique(z_aux), 
           cex = 0.5)
  }
  
  if (plt_univ) {
    readline(prompt="Plot, univariate: Press [enter] to continue")
    
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
      obs_univ <- obs_sub[, dd]
      plot.new()
      plot.window(ylim = c(min(obs_univ) - 0.3, max(obs_univ) + 0.3), 
                  xlim = c(1, N_sub),  xlab = "Time",  ylab = paste("y_", dd, sep=""), 
                  cex.lab = 1.1)
      #y_name = names(data_gesture)[1:D][dd]
      # y_name = rep("", D) # need to fix this
      title(xlab = "Time", ylab = paste("y_", dd,sep=""))
      for (j in 1:(length(cp_loc) - 1)) { 
        rect(cp_loc[j], ceiling(min(obs_univ)-2), cp_loc[j+1], ceiling(max(obs_univ) + 2), 
             col =scales::alpha(z_aux[cp_loc[j+1]], 0.2), border = F)
      }
      lines(1:N_sub, obs_univ, pch = 20, cex = 0.6, col = col_all[dd], type = "o")
      axis(1); axis(2)
      
      legend_aux <- numeric(K)
      for (kk in 1:K) {
        legend_aux[kk] <- paste(unique(state_seq_sub)[kk], 
                                unique(state_phase)[kk], sep = ": ")
      }
      legend("topright", 
             legend  = legend_aux, 
             pch = 15, 
             col = unique(z_aux), 
             cex = 0.5)
    }
    
  }
  
  return(list(obs = obs_sub, state = state_seq_sub, 
              state_phase = state_phase_sub))
}

Gesture_Phase_plot <- function(Gesture_Phase, only_acc = TRUE, plt_all = TRUE, plt_univ = FALSE)
{
  obs <- Gesture_Phase$obs
  D <- ncol(obs)
  
  if (only_acc) {obs <- obs[, (D-7):D]}
  
  state_seq <- Gesture_Phase$state
  state_phase <- Gesture_Phase$state_phase
  K <- length(unique(state_seq))
  
  #obs_sub <- obs[interval[1]:interval[2], ]
  #state_seq_sub <- state_seq[interval[1]:interval[2]]
  #state_phase_sub <- state_phase[interval[1]:interval[2]]
  #N_sub <- nrow(obs_sub)
  #D <- ncol(obs_sub)
  obs_sub <- obs
  state_seq_sub <- state_seq
  state_phase_sub <- state_phase
  N_sub <- nrow(obs_sub)
  D <- ncol(obs_sub)
  
  K <- length(unique(state_seq))
  if (K>2) {
    z_col <- brewer.pal(n = K, name = "Accent")
  } else {
    z_col <- c( "lightblue1", "lightsalmon")
  }
  
  z_aux <- sapply(1:N_sub, function(n) z_col[state_seq_sub[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N_sub)
  col_all <- rainbow(D)
  
  if (plt_all) {
    
    
    # - multivariate single plot
    par(mfrow = c(1, 1))
    title_main <- paste("Patient ", patient, "; Story ", story, sep = "")
    plot(1:N_sub, type = "n",  ylab =  "Sqrt(Physical Activity)", xlab = "Time",
         cex = 0.5, las = 1,
         xlim = c(1, N_sub),
         ylim = c(min(obs_sub) - 1, max(obs_sub) + 1), cex.lab = 1.1, 
         main = title_main)
    
    for (j in 1:(length(cp_loc) - 1)) { 
      rect(cp_loc[j], ceiling(min(obs_sub) - 2) , cp_loc[j+1], ceiling(max(obs_sub) + 2) , 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.4), border = F)
    }
    for (d in 1:D) {
      lines(obs_sub[, d], pch = 20,
            col = col_all[d], cex = 0.6, type = "o")
    }
    
    legend_aux <- numeric(K)
    for (kk in 1:K) {
      legend_aux[kk] <- paste(unique(state_seq_sub)[kk], 
                              unique(state_phase)[kk], sep = ": ")
    }
    legend("topright", 
           legend  = legend_aux, 
           pch = 15, 
           col = unique(z_aux), 
           cex = 0.5)
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
      obs_univ <- obs_sub[, dd]
      plot.new()
      plot.window(ylim = c(min(obs_univ) - 0.3, max(obs_univ) + 0.3), 
                  xlim = c(1, N_sub),  xlab = "Time",  ylab = paste("y_", dd, sep=""), 
                  cex.lab = 1.1)
      #y_name = names(data_gesture)[1:D][dd]
      # y_name = rep("", D) # need to fix this
      title(xlab = "Time", ylab = paste("y_", dd,sep=""))
      for (j in 1:(length(cp_loc) - 1)) { 
        rect(cp_loc[j], ceiling(min(obs_univ)-2), cp_loc[j+1], ceiling(max(obs_univ) + 2), 
             col =scales::alpha(z_aux[cp_loc[j+1]], 0.2), border = F)
      }
      lines(1:N_sub, obs_univ, pch = 20, cex = 0.6, col = col_all[dd], type = "o")
      axis(1); axis(2)
      
      legend_aux <- numeric(K)
      for (kk in 1:K) {
        legend_aux[kk] <- paste(unique(state_seq_sub)[kk], 
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



