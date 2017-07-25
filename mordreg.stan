data {
  int N;          //num obs
  int P;          //num predictors
  int D;          //num dimensions
  int L[D];       //num levels for each dimension
  
  int Y[N,D];  
  matrix[N,P] X;
}

parameters {
  matrix[P,D] beta;
  vector[D] alpha;
  vector<lower=0>[sum(L)-2*D] cutdiff;
}

transformed parameters {
  vector[sum(L)-D] cutpoint;
  
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
}

model {
  
  int start = 1;
  for (d in 1:D) {
    vector[L[d]-1] cpd;
    cpd = cutpoint[start:(start+L[d]-2)];
    start = start + L[d]-1;
    
    for (i in 1:N) Y[i,d] ~ ordered_logistic(alpha[d] + X[i,]*beta[,d],cpd);
  }
  
  to_vector(beta) ~ normal(0,5);
  alpha ~ normal(0,5);
  cutdiff ~ normal(0,5);
}

