# sparseVARHSMM
Sparse vector autoregressive (VAR) approximate hidden semi-Markov model (HSMM) 

**sparseVARHSMM** includes the stan software (and R utilities) to model temporal and contemporaneous (e.g. spatial) dependencies in multivariate time series data using a VAR HSMM, where the HSMM's generic state distribution is embedded in a special transition matrix structure (facilitating efficient likelihood evaluations and arbitrary approximation accuracy as in \citet{hadj2022bayesian}). Sparsity in the model is induced via the choice  the novel $l_1$-ball projection prior priors on the VAR coefficents. Such a prior introduces unconstrained latent variables and transforms them onto the space of VAR parameters in a way that provides positive probability that any element is exactly zero. Importantly,  
the transformation is almost surely continuous and differentiable, which allows for efficient posterior sampling using algorithms such as Hamiltonian Monte Carlo. We also include option for non-local priors on the parameters of the HSMM dwell distribution improving the ability of Bayesian model selection to distinguish whether the data is better supported by the simpler HMM, or more flexible HSMM. 

In the application of this research, we consider multivariate time series data that arise from a study on human gesture phase segmentation based on sensor data. As a segmentation exercise, I aim to model the data to identify periods of rest and active gesturing. 

