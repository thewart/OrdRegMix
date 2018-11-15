data {
  int N;          //num obs
  //int P;          //num predictors
  int D;          //num dimensions
  int K;          //num behavioral states

  int Y[N,D];
  real eta[N,D];

}

parameters {
  simplex[K] pi;		                  //topic weights
  matrix[D,K] gamma;                  //behavioral state effect
}

transformed parameters {
  vector[K] logp[N];

  for (i in 1:N) {
    logp[i] = log(pi);
    for (d in 1:D) {
        for (k in 1:K) logp[i,k] = logp[i,k] + bernoulli_lpmf(Y[i,d] | Phi(gamma[d,k] + eta[i,d]));
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
