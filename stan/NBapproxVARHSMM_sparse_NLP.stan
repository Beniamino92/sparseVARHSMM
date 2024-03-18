
functions {
  
   //Nonlocal Priors for High-Dimensional Estimation - David Rossella and Donatello Telesca
   //Tractable Bayesian Variable Selection: Beyond Normality - David Rossell & Francisco J. Rubio
   real pMOM_lpdf(real theta, real sigma2, real tau){
      return log((theta^2)) - log(tau*sigma2) + normal_lpdf(theta | 0, sqrt(tau*sigma2));
   }
   real eMOM_lpdf(real theta, real sigma2, real tau){
      return (sqrt(2) - tau*sigma2/theta^2) + normal_lpdf(theta | 0, sqrt(tau*sigma2));
   }
   // This is the heaviest tailed one
   real iMOM_lpdf(real theta, real sigma2, real tau){
      return 0.5*log(tau*sigma2)- 0.5*log(pi()) - 2*log(theta) - tau*sigma2/(theta^2);
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
  
  // phi = size of NB, phi = infty recovers the Poisson - Inverse gamma prior
  // phim1 = 1/size of NB, phim1 = 0 recovers the Poisson. - Gamma prior
  
     /** c_harzard dwell Negative-Binomial
  */
  real c_hazard_dwell_negativeBinomial(int r, real lambda, real phi){
    // shifted so that the minimum is 1
    return exp(neg_binomial_2_lpmf(r-1 | lambda, phi))/(1-neg_binomial_2_cdf(r-2, lambda, phi));
  }
  // E[NB] = lambda, var[NB] = lambda + lambda^2/phi
  
  /** B_ii Negative-Binomial - sparse representation of the transpose;
  */
  vector w_B_transpose_matrix_negativeBinomial_sparse(int K, int[] m_vect, matrix a_mat, vector lambda_vec, vector phi_vec)
  {
      int sum_m = sum(m_vect);
      vector[sum_m*K] w_B;// we can think about exactly getting the lengths right and indexing in if needs be
      int counter = 0;

      for(j in 1:K){
      // j_m is either 1, in which case we have the probabilities of arriving here from any other state 
        for(i in 1:K){
          if(i == j){
            // When m_j = 1 need to include the transition to itself - should this go in the id i == j?
            if(m_vect[j] == 1){
              counter += 1;
              w_B[counter] = (1 - c_hazard_dwell_negativeBinomial(m_vect[j], lambda_vec[j], phi_vec[j]));
            }
          } else{
            for(i_m in 1:m_vect[i]){
              counter += 1;
              w_B[counter] = (a_mat[i,j] * c_hazard_dwell_negativeBinomial(i_m, lambda_vec[i], phi_vec[i]));
            }
          }
        }
        if(m_vect[j] > 2){
        // j_m is in the middle, in which case we have it can only be reaced by the end state 
          for(j_m in 2:(m_vect[j] - 1)){
            counter += 1;
            w_B[counter] = (1 - c_hazard_dwell_negativeBinomial(j_m - 1, lambda_vec[j], phi_vec[j]));
          }
        }
        if(m_vect[j] > 1){
        // j_m = m_j and then it can be reached by the previous state or itself. 
          counter += 1;
          w_B[counter] =(1 - c_hazard_dwell_negativeBinomial(m_vect[j] - 1, lambda_vec[j], phi_vec[j]));
          counter += 1;
          w_B[counter] =(1 - c_hazard_dwell_negativeBinomial(m_vect[j], lambda_vec[j], phi_vec[j]));
        }
      }
      return w_B;
    }
    
  int[] v_B_transpose_matrix_negativeBinomial_sparse(int K, int[] m_vect, matrix a_mat, vector lambda_vec, vector phi_vec)
  {
      int sum_m = sum(m_vect);
      int v_B[K*sum_m];
      int counter = 0;
      
      int m_vect_temp[K+1];
      m_vect_temp[1] = 0;
      for(i in 1:K){
        m_vect_temp[i+1] = m_vect[i];//Adding a 0 to this makes the indexing below easier
      }
      
      for(j in 1:K){
      // j_m is either 1, in which case we have the probabilities of arriving here from any other state 
        for(i in 1:K){
          if(i == j){
            // When m_j = 1 need to include the transition to itself - should this go in the id i == j?
            if(m_vect_temp[j + 1] == 1){
              counter += 1;
              v_B[counter] = (sum(m_vect_temp[1:j]) + m_vect_temp[j + 1]);
            }
          } else{
            for(i_m in 1:m_vect_temp[i + 1]){
              counter += 1;
              v_B[counter] = (sum(m_vect_temp[1:i]) + i_m);
            }
          }
        }
        if(m_vect_temp[j + 1] > 2){
        // j_m is in the middle, in which case we have it can only be reaced by the end state 
          for(j_m in 2:(m_vect_temp[j + 1] - 1)){
            counter += 1;
            v_B[counter] = (sum(m_vect_temp[1:j]) + j_m - 1); 
          }
        }
        if(m_vect_temp[j + 1] > 1){
        // j_m = m_j and then it can be reached by the previous state or itself. 
          counter += 1;
          v_B[counter] =(sum(m_vect_temp[1:j]) + m_vect_temp[j + 1] - 1);
          counter += 1;
          v_B[counter] =(sum(m_vect_temp[1:j]) + m_vect_temp[j + 1]);
        }
      }
      return v_B;
    
    }
  
  int[] u_B_transpose_matrix_negativeBinomial_sparse(int K, int[] m_vect, matrix a_mat, vector lambda_vec, vector phi_vec)
  {
      int sum_m = sum(m_vect);
      int u_B[sum_m + 1];
      int counter = 0;
      int ind;
      
      counter += 1;
      u_B[counter] = 1;
      for(j in 1:K){
      // j_m is either 1, in which case we have the probabilities of arriving here from any other state 
        ind = 0;
        for(i in 1:K){
          if(i == j){
            // When m_j = 1 need to include the transition to itself - should this go in the id i == j?
            if(m_vect[j] == 1){
              ind += 1;
            }
          } else{
            for(i_m in 1:m_vect[i]){
              ind += 1;
            }
          }
        }
        counter += 1;
        u_B[counter] = (u_B[counter - 1] + ind);
      if(m_vect[j] > 2){
      // j_m is in the middle, in which case we have it can only be reaced by the end state 
        for(j_m in 2:(m_vect[j] - 1)){
          counter += 1;
          u_B[counter] = (u_B[counter - 1] + 1); // c[1:size(c)]
        }
      }
      if(m_vect[j] > 1){
        // j_m = m_j and then it can be reached by the previous state or itself. 
        counter += 1;
        u_B[counter] = (u_B[counter - 1] + 2);
      }
    }
      return u_B;
  }
    
  /** convert (K-1) simplex gamma to a K*K matrix transition matrix.
  */
  matrix a_to_mat(int K, vector[] a)  
  {
    matrix[K, K] A = rep_matrix(0, K, K); 
    int count_j = 0;
    for (i in 1:K) {
      count_j = 0;
      for (j in 1:K) {
        if(i != j) {
          count_j += 1;
          A[i, j] = a[i][count_j]; 
        }
      }
    }
    return A;
  }
  
  
  /**  generate a N * M matrix with
  multivariate Gaussian emissions for t=1,..,N and j=1,..,M=sum(m).
  */
  matrix log_VAR_Emissions(int N, int K, int D, int P, int[] m, vector[] y, vector[] beta_0, 
                       vector[] beta, vector[] tau, matrix[] L_Omega) 
  {
    int M = sum(m);
    matrix[N, M] all_logprobs;
    vector[D] mu;
    matrix[D, D] beta_temp;
    
    matrix[D, D] Omega;
    matrix[D, D] Sigma[K];
    for(k in 1:K){
      Omega = multiply_lower_tri_self_transpose(L_Omega[k]);
      Sigma[k] = quad_form_diag(Omega, tau[k]);
    }
    
    
    // first P terms = 1
    for (k in 1:M) {
      for (n in 1:P) {
              //allprobs[n, k] = exp(multi_normal_lpdf(y[n] | y[n], 
              //                     Sigma[k]));
              all_logprobs[n, k] = 0;
                                  
      }
    }
    
    // VAR emissions: (P+1):N, 
    
    for (n in (P+1):N) {
      for(k in 1:K) {
        mu = beta_0[k];
        for (p in 1:P) {
          beta_temp = get_beta_temp(beta[k], p, D);
          mu += beta_temp * y[n - p];
        }
        if(k ==1) {
          for (i in 1:m[1]) {
            all_logprobs[n, i] = multi_normal_lpdf(y[n] | mu, 
                             Sigma[1]);
          }
        }
        else {
          for(i in 1:m[k]) {
            all_logprobs[n, sum(m[1:(k-1)]) + i] = multi_normal_lpdf(y[n] | mu, 
                                                   Sigma[k]);
          }
        }
      }
    }
  
    return all_logprobs; 
  }
  
  
  /** perform forward algorithm to 
  generate alpha_t (j) for t=1,..,N and j=1,.., M = sum(m)
  via forward dynamic programming - similar to p.243 (Zucchini et al.).
  */
  real forward_messages(int N, int M, matrix log_emissions, vector w_B_transpose, int[] v_B_transpose, int[] u_B_transpose) 
  {
    vector[M] log_foo;
    vector[M] log_foo_temp;
    real log_sumfoo;
    real lscale;
    // alpha_1
    for (k in 1:M) {
      log_foo_temp[k] = log_emissions[1, k];
    }
    log_sumfoo = log_sum_exp(log_foo_temp);
    lscale = log_sumfoo;
    log_foo = log_foo_temp - log_sumfoo;
    
    // alpha_t, t = 2, ..., N
    for (i in 2:N) {
      for(k in 1:M){

        if(u_B_transpose[k+1] == u_B_transpose[k]){
          log_foo_temp[k] = negative_infinity();
        } else{
          log_foo_temp[k] = log_sum_exp(log_foo[v_B_transpose[u_B_transpose[k]:(u_B_transpose[k+1]-1)]] + log(w_B_transpose[u_B_transpose[k]:(u_B_transpose[k+1]-1)])) + log_emissions[i, k];// this works amazingly, faster and same inference!
        }
      }
      log_sumfoo = log_sum_exp(log_foo_temp);
      lscale += log_sumfoo;
      log_foo = log_foo_temp - log_sumfoo;
    }
    return lscale;
  }

  
  /** compute likelihood p(y_1, .., y_T  | )
  */
  real llk_lp(int N, int K, int D, int P, int[] m, vector[] y, vector[] beta_0,
              vector[] beta, vector[] tau, matrix[] L_Omega, 
              vector lambda, vector phi, vector[] a)
  {
    int M = sum(m);
    matrix[K, K] A =  a_to_mat(K, a);
    vector[K*M] w_B_transpose = w_B_transpose_matrix_negativeBinomial_sparse(K, m, A, lambda, phi);
    int v_B_transpose[K*M] = v_B_transpose_matrix_negativeBinomial_sparse(K, m, A, lambda, phi);
    int u_B_transpose[M + 1] = u_B_transpose_matrix_negativeBinomial_sparse(K, m, A, lambda, phi);
    matrix[N, M] log_allprobs = log_VAR_Emissions(N, K, D, P, m, y, beta_0, beta, tau, L_Omega);
    real llk =  forward_messages(N, M, log_allprobs, w_B_transpose, v_B_transpose, u_B_transpose);
    return llk;
  } 
}

data {
  int<lower=0> N; // length time series
  int<lower=1> D; // dimensionality // lower bound is 1 becuase of the deomposition of mu
  int<lower=1> P; // lag of the VAR
  int<lower=0> K; // number of states
  int m[K]; // size state aggr
  vector[D] y[N]; // data

  
  // hyperparms
  vector[D] beta_0_loc[K]; 
  real<lower=0> beta_0_scale; 
  real beta_loc; 
  real<lower=0> beta_scale;
  real tau_loc; 
  real<lower=0> tau_scale;
  real<lower=0> Omega_shape;
  
  real<lower = 0> a_0_lambda[K];  // prior rate NB dwell
  real<lower = 0> b_0_lambda[K]; //         ""
  real<lower = 0> tau_phi[K];  //           ""
  real<lower = 0> sigma2_phi[K]; //            ""
  
  vector<lower = 0>[K-1] alpha_0[K]; // prior dirichlet probs
  
}

parameters {
  simplex[K-1] gamma[K];// transition prob mat 
  vector<lower=0, upper=1>[K] tilde_p; // lambda = p/(1-p)
  vector<lower=0, upper=1>[K] tilde_phim1; // phim1 = tilde_phim1/(1-tilde_phim1)
  vector[D] beta_0[K]; // mean paramaters
  vector[D*D*P] beta[K]; // AR parameters
  vector<lower=0>[D] tau[K];  // scales gauss emission
  cholesky_factor_corr[D] L_Omega[K]; // cholesky factor correlation
}


model {
  
  vector[K] phi;
  vector[K] lambda;
  vector[K] phim1;
  
  for(k in 1:K){
    lambda[k] = tilde_p[k]/(1-tilde_p[k]);
    phim1[k] = tilde_phim1[k]/(1-tilde_phim1[k]);
    // transforming phi
    phi[k] = 1/phim1[k];
  }

  // priors
  for(k in 1:K){
    target += gamma_lpdf(tilde_p[k]/(1-tilde_p[k]) | a_0_lambda, b_0_lambda) - 2*log(1-tilde_p[k]);
    target += eMOM_lpdf(log(tilde_phim1[k]/(1-tilde_phim1[k])) | sigma2_phi[k], tau_phi[k]) + log(1/tilde_phim1[k] + 1/(1-tilde_phim1[k]));
    target += lkj_corr_cholesky_lpdf(L_Omega[k] | Omega_shape);
    target += cauchy_lpdf(tau[k] | tau_loc, tau_scale);
    target += dirichlet_lpdf(gamma[k] | alpha_0[k]);
    target += double_exponential_lpdf(beta_0[k] | beta_0_loc[k], beta_0_scale);
    target += double_exponential_lpdf(beta[k] | beta_loc, beta_scale);
  }
  // likelihood
  target += llk_lp(N, K, D, P, m, y, beta_0, beta, tau,  
                   L_Omega, lambda, phi, gamma); 
}


generated quantities {
  vector[K] lambda;
  vector[K] phim1;
  
  for(k in 1:K){
    lambda[k] = tilde_p[k]/(1-tilde_p[k]);
    phim1[k] = tilde_phim1[k]/(1-tilde_phim1[k]);
  }
  
}




