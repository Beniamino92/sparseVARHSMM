# sparseVARHSMM
###  <ins>Sparse vector autoregressive (VAR) approximate hidden semi-Markov model (HSMM)<ins>

**sparseVARHSMM** includes the *stan* software (and R utilities) to model temporal and contemporaneous (e.g. spatial) dependencies in multivariate time series data using a VAR HSMM, where the HSMM's generic state distribution (e.g. negative-binomial) is embedded in a ***special transition matrix*** structure (facilitating efficient likelihood evaluations and arbitrary approximation accuracy, as in [BayesApproxHSMM](https://github.com/Beniamino92/BayesianApproxHSMM/)). ***Sparsity*** in the model is induced via the choice  the novel $l_1$-ball projection prior priors on the VAR coefficents. Such a prior introduces unconstrained latent variables and transforms them onto the space of VAR parameters in a way that provides positive probability that any element is exactly zero. Importantly, the transformation is almost surely continuous and differentiable, which allows for efficient posterior sampling using algorithms such as Hamiltonian Monte Carlo. We also include options for ***non-local priors*** (NLP) on the parameters of the HSMM dwell distribution improving the ability of Bayesian model selection to distinguish whether the data is better supported by the simpler HMM, or more flexible HSMM. 

This software allows for the following four modelling option 

1. VAR HMM
2. VAR HSMM [ $l_1$-ball ]
3. VAR HSMM [ NLP ]
4. VAR HSMM [ $l_1$-ball  and NLP ]
  
You can find below some representative examples of our software facilities, applied to real sensor data, where we consider an approximate VAR HSMM [ $l_1$-ball  and NLP ] with negative-binomial dwell durations. Different dwell/emission distributions from the ones considered in this paper can be easily developed. Users need only change the corresponding function in our stan files (and R utilities). 
  

```r
NBapproxVARHSMM_stan <- stan_model(file = "stan/NBapproxVARHSMM_sparse_l1ball_fullCov_NLP.stan")
fit <- sampling(object = NBapproxVARHSMM_stan,
                data = NBapproxVARHSMM_data, seed = 123, 
                chains = 1, iter = 1000 + N_MCMC, 
                warmup = 1000)  
```

  say about viterbi algorithm. 
```r
predictive <- NBapproxVARHSMM_getPredictive(fit , m, obs, pseudo = FALSE, 
                                            L1_ball = TRUE, ndraw = 50)
z_hat <- predictive$z_hat
y_hat <- predictive$y_hat
plotPosteriorPredictive(obs, y_hat, z_hat, K)
```
         
<p align="center">
<img src="https://github.com/Beniamino92/sparseVARHSMM/blob/main/figures/postpred_training.png" width="700" heigth="100"/> 
</p>
  
```r
 # state probabilities plots 
 NBapproxVARHSMM_plotStateProbs(fit, obs, m = m, pseudo = FALSE, L1_ball = TRUE)
```
  
<p align="center">
<img src="https://github.com/Beniamino92/sparseVARHSMM/blob/main/figures/stateprobs_training.png" width="700" heigth="100"/> 
</p>
  
  
```r
params <- rstan::extract(fit)
p_est <- get_inclusion_sims(params, p_est = TRUE)$p_est
beta_est <- get_beta_est(params, L1_ball = TRUE, mat = TRUE)
```
  
```r
 # Rest
plotDAG(p_est[1, , , 1], ylabels, color = "lightblue1", main = "Rest")
# Active
plotDAG(p_est[2, , , 1], ylabels, color = "lightsalmon", main = "Active")
```
  
<p align="center">
<img src="https://github.com/Beniamino92/sparseVARHSMM/blob/main/figures/DAGactive.png" width="600" heigth="600"/> 
</p>
  
```r
 # PPI (Rest)
plot.heat(Matrix = p_est[1, , , 1], Xlab="", Ylab="", Main="PPI (Rest)", limit=c(0,1))
# VAR coeffs (Active)
plot.heat(Matrix = beta_est[2, , , 1], Xlab="", Ylab="", Main="VAR coeffs (Active)", limit=c(-1,1))
```

<p align="center">
<img src="https://github.com/Beniamino92/sparseVARHSMM/blob/main/figures/inclusion_coeffs_covariance.png" width="600" heigth="600"/> 
</p>
  

<!-- In the application of this research, we consider multivariate time series data that arise from a study on human gesture phase segmentation based on sensor data. As a segmentation exercise, We aim to model the data to identify periods of rest and active gesturing.  -->

