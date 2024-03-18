functions {
  

   vector l1_ball_projection(vector beta, real r) {
     if(sum(fabs(beta)) <= r){return beta;}
     else{
      int p = num_elements(beta);
      vector[p] abs_beta_sort = sort_desc(fabs(beta));
      vector[p] cumsum_abs_beta_sort = cumulative_sum(abs_beta_sort);
      vector[p] mu;
      vector[p] theta;
      //int c = 0;
      int c = p;
      for(j in 1:p){
         mu[j] = max({cumsum_abs_beta_sort[j] - r, 0});
         if(abs_beta_sort[j] > mu[j]/j){
            c = j;
         }
      }
      for(i in 1:p){
         theta[i] = fabs(beta[i])/beta[i]*max({fabs(beta[i]) - mu[c]/c, 0});
      }      
      return theta;
     }
    }
  
  
  /**
  !
  */
  matrix get_beta_temp(vector beta, int p, int D) 
  {
    int start;
    int end; 
    if (p == 1) {
      start = 1;
    } else {
      start = ((D*D)*(p-1)) + 1;
    }
    end = (D*D)*p;
    return(to_matrix(beta[start:end], D, D, 0));
  }
  
  /**  generate a N * K matrix with
  multivariate Gaussian emissions for t=1,..,N and j=1,..,K.
  */
  matrix VAR_Emissions(int N, int K, int D, int P, vector[] y, vector[] beta_0, 
                       vector[] beta, vector[] tau, matrix[] L_Omega) 
  {
    //matrix[N, K] allprobs;
    matrix[N, K] alllogprobs;
    vector[D] mu;
    matrix[D, D] beta_temp;
    
    matrix[D, D] Omega;
    matrix[D, D] Sigma[K];
    for(k in 1:K){
      Omega = multiply_lower_tri_self_transpose(L_Omega[k]);
      Sigma[k] = quad_form_diag(Omega, tau[k]);
    }
    

    for (k in 1:K) {
      for (n in 1:P) {
              alllogprobs[n, k] = 0;
                                  
      }
    }
  
    
    for (k in 1:K) {
      for (n in (P+1):N) {
        mu = beta_0[k];
        for (p in 1:P) {
          beta_temp = get_beta_temp(beta[k], p, D);
          mu += beta_temp * y[n - p];
        }

        alllogprobs[n, k] = multi_normal_lpdf(y[n] | mu, 
                             Sigma[k]);
      }
    }
    //return allprobs; 
    return alllogprobs;
  }
  
    /** perform forward algorithm to 
  generate alpha_t (j) for t=1,..,N and j=1,.., M = sum(m)
  via forward dynamic programming - p.243 (Zucchini et al.).
  */
  //real forwardMessages(int N, int K, matrix emissions, matrix gamma_mat)
  real forwardMessages(int N, int K, matrix log_emissions, matrix gamma_mat)
  {
    //vector[K] foo;
    vector[K] log_foo;
    //real sumfoo;
    real lscale;
    vector[K] gamma_init; 
    // alpha_1
    for (k in 1:K) {
      log_foo[k] = log_emissions[1, k];
    }
    lscale = log_sum_exp(log_foo);
    log_foo = log_foo - lscale;
    // alpha_t, t = 2, ..., N
    for (i in 2:N) {
      log_foo = (log(exp(log_foo)'*gamma_mat) + log_emissions[i, :])';
      lscale = lscale + log_sum_exp(log_foo);
      log_foo = log_foo - log_sum_exp(log_foo);
    }
    return lscale;
  }
  
  /** convert simplex gamma to a K*K matrix. 
  */
  matrix gammaToMat(int K, vector[] gamma)
  {
    matrix[K, K] gamma_mat;
    // Creating Gamma Matrix for Forward
    for(i in 1:K)
      for(j in 1:K)
        gamma_mat[i, j] = gamma[i][j];
    return gamma_mat;
  }

  /** compute likelihood p(y_1, .., y_T  | )
  */
  real llk_lp(int N, int K, int D, int P, vector[] y, vector[] beta_0,
              vector[] beta, vector[] tau, matrix[] L_Omega,  vector[] gamma)
  {
    matrix[N, K] emissions = VAR_Emissions(N, K, D, P, y, beta_0, beta, tau, L_Omega);
    matrix[K, K] gamma_mat = gammaToMat(K, gamma);
    real llk =  forwardMessages(N, K, emissions, gamma_mat);
    return llk;
  }
}

data {
  int<lower=1> N; // length time series
  int<lower=1> D; // dimensionality // lower bound is 1 becuase of the deomposition of mu
  int<lower=1> K; // number of states
  int<lower=1> P; // lag of the VAR
  vector[D] y[N]; // data

  // hyperparms
  vector[D] beta_0_loc[K]; 
  real<lower=0> beta_0_scale; 
  real beta_loc; 
  real<lower=0> beta_scale;
  real tau_loc; 
  real<lower=0> tau_scale;
  real<lower=0> Omega_shape;
  vector<lower = 0>[K-1] alpha_0_rest[K]; // prior dirichlet probs
  real<lower=0> alpha_0_diag[K];
  real<lower=0> r_alpha;
}

parameters {
  vector<lower=0, upper=1>[K] gamma_diag;
  simplex[K-1] gamma_rest[K];// transition prob mat 
  vector[D] beta_0[K]; // mean gauss emission - rest of the dimension 
  vector[D*D*P] tilde_beta[K]; // these are scaled by $N$ to get convergence
  vector<lower=0>[D] tau[K];  // scales gauss emission
  cholesky_factor_corr[D] L_Omega[K]; // additional if we want to model conditional covariance 
  vector<lower=0>[K] tilde_ray;
}


model {
  
  vector[K] ray;
  vector[D*D*P] beta[K];
  vector[D*D*P] theta[K];
  vector[K] gamma[K];
  for (k in 1:K) {
    ray[k] = tilde_ray[k]*sqrt(D*D*P); // we cancle the denominator cause r is the sume of all of these 
    beta[k] = tilde_beta[k]*sqrt(N);
    theta[k] = l1_ball_projection(beta[k], ray[k]);
  }
  for(k in 1:K){
    for(k_prime in 1:K){
      if(k == k_prime){
        gamma[k, k] = gamma_diag[k];
      } 
      if(k > k_prime){
        gamma[k, k_prime] = gamma_rest[k, k_prime]*(1-gamma_diag[k]);
      }
      if(k < k_prime){
        gamma[k, k_prime] = gamma_rest[k, k_prime - 1]*(1-gamma_diag[k]);
      }
      
    }
  }
  
  // priors
  for (k in 1:K) {
    target += lkj_corr_cholesky_lpdf(L_Omega[k] | Omega_shape);
    target += cauchy_lpdf(tau[k] | tau_loc, tau_scale);
    target += beta_lpdf(gamma_diag[k] | alpha_0_diag[k], sum(alpha_0_rest[k]));
    target += dirichlet_lpdf(gamma_rest[k] | alpha_0_rest[k]);
    target += double_exponential_lpdf(beta_0[k] | beta_0_loc[k], beta_0_scale);
    target += double_exponential_lpdf(tilde_beta[k] | beta_loc, beta_scale/sqrt(N));
    target += exponential_lpdf(tilde_ray[k] | r_alpha*sqrt(D*D*P));

  }

  target += llk_lp(N, K, D, P, y, beta_0, theta, tau, L_Omega, gamma); 
}



//
generated quantities {
  vector[K] ray;
  vector[D*D*P] beta[K];
  vector[D*D*P] theta[K];
  matrix[K, K] gamma;
  for (k in 1:K) {
    ray[k] = tilde_ray[k]*sqrt(D*D*P);
    beta[k] = tilde_beta[k]*sqrt(N);
    theta[k] = l1_ball_projection(beta[k], ray[k]);
  }
  for(k in 1:K){
    for(k_prime in 1:K){
      if(k == k_prime){
        gamma[k, k] = gamma_diag[k];
      } 
      if(k > k_prime){
        gamma[k, k_prime] = gamma_rest[k, k_prime]*(1-gamma_diag[k]);
      }
      if(k < k_prime){
        gamma[k, k_prime] = gamma_rest[k, k_prime - 1]*(1-gamma_diag[k]);
      }
      
    }
  }
}


