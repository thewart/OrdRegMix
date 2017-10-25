functions {
    vector ordered_logistic_probabilities(real eta, vector c) {
      int L = num_elements(c) + 1;
      vector[L] theta;

      theta[1] = inv_logit(c[1]-eta);
      for (l in 2:(L-1)) theta[l] = inv_logit(c[l]-eta) - inv_logit(c[l-1]-eta);
      theta[L] = 1 - inv_logit(c[L-1]-eta);
      return theta;
    }
    
    real ordered_logistic_grouped_lpmf(int[] y, real eta, vector c) {
      int L = num_elements(c) + 1;
      vector[L] theta = ordered_logistic_probabilities(eta,c);
      return multinomial_lpmf(y | theta);
    }
    
    real ordered_logistic_expectation(real eta, vector c) {
      real ex = 0.0;
      int L = num_elements(c) + 1;
      vector[L] theta = ordered_logistic_probabilities(eta,c);
      
      for (l in 1:L) ex = ex + theta[l]*l;
      return ex;
    }
}

data {
  int N;          //num obs
  int P1;         //num predictors
  int P2;         //num regularized predictors
  int D;          //num dimensions
  int L[D];       //num levels for each dimension
  int V;          //number of variance components
  int K[V];       //number of levels per variance component
  
  int Y[N,sum(L)];
  matrix[N,P1] X1;
  matrix[N,P2] X2;
  int Z[N,V];
  
  real<lower=1> nu; //degrees of freedom for gamma
}

transformed data {
  int vstart[V+1];

  vstart[1] = 1;
  vstart[V+1] = sum(K)+1;
  for (v in 2:V) vstart[v] = vstart[v-1] + K[v-1];
}

parameters {
  real<lower=0> sigma_u[V,D];
  real<lower=0> sigma_l;
  vector[D] alpha;                   //intercept
  matrix[P1,D] beta;                 //reg coeffs 
  matrix[sum(K),D] u_raw;            //raw random effects
  matrix[P2,D] lambda_raw;           //raw regularized reg coeffs
  vector<lower=0>[sum(L)-2*D] cutdiff;
}

transformed parameters {
  vector[sum(L)-D] cutpoint;
  vector[N] eta[D];
  matrix[P2,D] lambda = sigma_l*lambda_raw;
  
  {
    int start = 1;

    for (d in 1:D) {
      
      cutpoint[start]=0;
      for (l in 1:(L[d]-2) ) {
        cutpoint[start+l] = cutpoint[start+l-1] + cutdiff[start+l-d];
      }
      start = start + L[d]-1;
      
      eta[d] = alpha[d] + X1*beta[,d] + X2*lambda[,d]; //+ sigma_g * L_A * g_raw;
      for (v in 1:V) {
        vector[K[v]] u;
        u = sigma_u[v,d] * u_raw[vstart[v]:(vstart[v+1]-1),d];
        for (i in 1:N) eta[d][i] = eta[d][i] + u[Z[i,v]];
      }
    }
  }
}

model {
  
  int start = 1;
  for (d in 1:D) {
    vector[L[d]-1] cpd;

    cpd = cutpoint[start:(start+L[d]-2)];
    for (i in 1:N) Y[i,(start+d-1):(start+d+L[d]-2)] ~ ordered_logistic_grouped(eta[d][i],cpd);
    start = start + L[d]-1;
  }
  alpha ~ normal(0,5);
  to_vector(beta) ~ normal(0,5);
  to_vector(u_raw) ~ normal(0,1);
  to_vector(lambda_raw) ~ student_t(nu,0,1);
  for (d in 1:D) sigma_u[,d] ~ normal(0,5);
  sigma_l ~ normal(0,5);
  cutdiff ~ normal(0,5);
}

generated quantities {
  real log_lik[N];
  real pseudo_r2[D];
  real r2[D];
  
  for (i in 1:N) log_lik[i] = 0;
  {
    int start = 1;
    for (d in 1:D) {
      vector[L[d]-1] cpd;
      real yhat_full;
      real yhat_red;
      real yhat_alpha;
      real rss_full;
      real rss_red;
      real rss_alpha;

      rss_full = 0;
      rss_red = 0;
      rss_alpha = 0;
      cpd = cutpoint[start:(start+L[d]-2)];
      for (i in 1:N) {
        #calc loglik
        log_lik[i] = log_lik[i] + 
        ordered_logistic_grouped_lpmf(Y[i,(start+d-1):(start+d+L[d]-2)] | eta[d][i], cpd);
        
        #calc rhats
        yhat_full = ordered_logistic_expectation(eta[d][i],cpd);
        yhat_red = ordered_logistic_expectation(eta[d][i]-X2[i,]*lambda[,d], cpd);
        yhat_alpha = ordered_logistic_expectation(alpha[d],cpd);
        for (l in 1:L[d]) {
          rss_full = rss_full + pow(yhat_full-l,2) * Y[i,(start+d-1):(start+d+L[d]-2)][l];
          rss_red = rss_red + pow(yhat_red-l,2) * Y[i,(start+d-1):(start+d+L[d]-2)][l];
          rss_alpha = rss_alpha + pow(yhat_alpha-l,2) * Y[i,(start+d-1):(start+d+L[d]-2)][l];
        }
      }
      pseudo_r2[d] = 1 - rss_full/rss_red;
      r2[d] = 1 - rss_full/rss_alpha;
      start = start + L[d]-1;
    }
  }
  
}
