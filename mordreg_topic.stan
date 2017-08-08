data {
  int N;          //num obs
  int P;          //num predictors
  int D;          //num dimensions
  int L[D];       //num levels for each dimension
  int K;          //num behavioral states
  //int V;          //number of variance components
  //int K[V];       //number of levels per variance component
  
  int Y[N,D];
  matrix[N,P] X;
  //int Z[N,V];
  
  //real<lower=1> nu; //degrees of freedom for gamma
}

// transformed data {
//   int vstart[V+1];
// 
//   vstart[1] = 1;
//   vstart[V+1] = sum(K)+1;
//   for (v in 2:V) vstart[v] = vstart[v-1] + K[v-1];
// }

parameters {
  simplex[K] pi;		                   //topic weights
  //real<lower=0> sigma_u[V,D];
  //vector[sum(K)] u_raw[D];           //raw random effects
  matrix[P,D] beta;                 //reg coeffs 
  //vector[D] alpha;                     //intercept
  matrix[D,K] gamma;                   //behavioral state effect
  vector<lower=0>[sum(L)-2*D] cutdiff;
}

transformed parameters {
  vector[sum(L)-D] cutpoint;
  //vector[N] eta[D];
  {
    int start = 1;

    for (d in 1:D) {
      
      cutpoint[start]=0;
      for (l in 1:(L[d]-2) ) {
        cutpoint[start+l] = cutpoint[start+l-1] + cutdiff[start+l-d];
      }
      start = start + L[d]-1;
      
      // eta[d] = alpha[d] + X1*beta[,d] + X2*lambda[,d]; //+ sigma_g * L_A * g_raw;
      // for (v in 1:V) {
      //   vector[K[v]] u;
      //   u = sigma_u[v,d] * u_raw[d][vstart[v]:(vstart[v+1]-1)];
      //   for (i in 1:N) eta[d][i] = eta[d][i] + u[Z[i,v]];
      // }
    }
  }
}

model {
  vector[K] logp[N];
  int start = 1;
  
  for (i in 1:N) logp[i] = log(pi);
  for (d in 1:D) {
    vector[L[d]-1] cpd;

    cpd = cutpoint[start:(start+L[d]-2)];
    for (i in 1:N) {
      for (k in 1:K) logp[i,k] = logp[i,k] + ordered_logistic_lpmf(Y[i,d] | gamma[d,k] + X[i,]*beta[,d],cpd);
      //ordered_logistic_lpmf(Y[i,d] | alpha[d] + gamma[d,k],cpd);
    }
    start = start + L[d]-1;
    
    // sigma_u[,d] ~ normal(0,2);
    // u_raw[d] ~ normal(0,1);
  }
  for (i in 1:N) target += log_sum_exp(logp[i]);
  
  //to_vector(beta) ~ student_t(4,0,2.5);
  //to_vector(gamma) ~ normal(0,1);
  //sigma_g ~ normal(0,2.5);
  //alpha ~ normal(0,2.5);
  //cutdiff ~ normal(0,5);
}

