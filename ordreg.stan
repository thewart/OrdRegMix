functions {
    real ordered_logistic_grouped_lpmf(int[] y, real eta, vector c) {
      int L = num_elements(c) + 1;
      vector[L] theta;

      theta[1] = inv_logit(c[1]-eta);
      for (l in 2:(L-1)) theta[l] = inv_logit(c[l]-eta) - inv_logit(c[l-1]-eta);
      theta[L] = 1 - inv_logit(c[L-1]-eta);
      return multinomial_lpmf(y | theta);
    }
}

data {
  int N;          //num obs
  int P;          //num predictors
  int L;          //num levels
  int V;          //number of variance components
  int K[V];       //number of levels per variance component
  //int R;          //number of biological replicates

  int y[N,L];  
  matrix[N,P] X;
  int Z[N,V];     //levels
  //matrix[N,R] L_A;//cholesky-like kinship matrix
}

transformed data {
  int vstart[V+1];
   
  vstart[1] = 1;
  vstart[V+1] = sum(K)+1;
  for (v in 2:V) vstart[v] = vstart[v-1] + K[v-1];
}

parameters {
  real<lower=0> sigma_u[V];
  //real<lower=0> sigma_g;
  vector[sum(K)] u_raw;           //raw random effects
  //vector[R] g_raw;                //raw genetic effects
  vector[P] beta;
  real alpha;
  vector<lower=0>[L-2] cutdiff;
}

transformed parameters {
  vector<lower=0>[L-1] cutpoint;
  vector[N] eta;
  
  cutpoint[1]=0;
  for (l in 1:(L-2) ) {
    cutpoint[l+1] = cutpoint[l] + cutdiff[l];
  }
  

  eta = alpha + X*beta; //+ sigma_g * L_A * g_raw;
  for (v in 1:V) {
    vector[K[v]] u;
    u = sigma_u[v] * u_raw[vstart[v]:(vstart[v+1]-1)];
    for (i in 1:N) eta[i] = eta[i] + u[Z[i,v]];
  }
}

model {
  
  for (i in 1:N) y[i] ~ ordered_logistic_grouped(eta[i],cutpoint);
  sigma_u ~ normal(0,2);
  u_raw ~ normal(0,1);
  //sigma_g ~ normal(0,2);
  //g_raw ~ normal(0,1);
  beta ~ normal(0,2);
  alpha ~ normal(0,2);
  cutdiff ~ normal(0,5);
}

