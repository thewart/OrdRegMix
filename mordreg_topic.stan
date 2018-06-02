functions {
    vector ordered_probit_probabilities(real eta, vector c) {
      int L = num_elements(c) + 1;
      vector[L] theta;

      theta[1] = Phi(c[1]-eta);
      for (l in 2:(L-1)) theta[l] = Phi(c[l]-eta) - Phi(c[l-1]-eta);
      theta[L] = 1 - Phi(c[L-1]-eta);
      return theta;
    }
    
    real ordered_probit_lpmf(int y, real eta, vector c) {
      int L = num_elements(c) + 1;
      vector[L] theta = ordered_probit_probabilities(eta,c);
      return categorical_lpmf(y | theta);
    }
}

data {
  int N;          //num obs
  //int P;          //num predictors
  int D;          //num dimensions
  int L[D];       //num levels for each dimension
  int K;          //num behavioral states

  int Y[N,D];
  real eta[N,D];

}

parameters {
  simplex[K] pi;		                  //topic weights
  matrix[D,K] gamma;                  //behavioral state effect
  vector<lower=0>[sum(L)-2*D] cutdiff;
}

transformed parameters {
  vector[sum(L)-D] cutpoint;
  vector[K] logp[N];

  {
    int start = 1;

    for (d in 1:D) {
      
      cutpoint[start]=0;
      for (l in 1:(L[d]-2) ) {
        cutpoint[start+l] = cutpoint[start+l-1] + cutdiff[start+l-d];
      }
      start = start + L[d]-1;
      
    }
  }
  
  {
    int start = 1;
  
    for (i in 1:N) logp[i] = log(pi);
    for (d in 1:D) {
      vector[L[d]-1] cpd;

      cpd = cutpoint[start:(start+L[d]-2)];
      for (i in 1:N)
        for (k in 1:K) logp[i,k] = logp[i,k] + ordered_probit_lpmf(Y[i,d] | gamma[d,k] + eta[i,d],cpd);

      start = start + L[d]-1;
    } 
  }
  
  
}

model {
  for (i in 1:N) target += log_sum_exp(logp[i]);
  
  //to_vector(gamma) ~ normal(0,1);
  //cutdiff ~ normal(0,5);
}

generated quantities {
  real r[N,K];
  
  for (i in 1:N) {
    real norm;
    norm = log_sum_exp(logp[i]);
    for (k in 1:K) r[i,k] = exp(logp[i,k] - norm);
  }
}
