data {
  int N;          //num obs
  int P;          //num predictors
  int L;          //num levels
  
  int y[N];  
  matrix[N,P] X;
}

parameters {
  vector[P] beta;
  real alpha;
  vector<lower=0>[L-2] cutdiff;
}

transformed parameters {
  vector<lower=0>[L-1] cutpoint;
  
  cutpoint[1]=0;
  for (l in 1:(L-2) ) {
    cutpoint[l+1] = cutpoint[l] + cutdiff[l];
  }
}

model {
  
  for (i in 1:N) y[i] ~ ordered_logistic(alpha + X[i,]*beta,cutpoint);
  
  beta ~ normal(0,5);
  alpha ~ normal(0,5);
  cutdiff ~ normal(0,5);
}

