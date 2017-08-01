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
  int D;          //num dimensions
  int L[D];       //num levels for each dimension
  int V;          //number of variance components
  int K[V];       //number of levels per variance component
  
  int Y[N,sum(L)];
  matrix[N,P] X;
  int Z[N,V];
}

transformed data {
  int vstart[V+1];

  vstart[1] = 1;
  vstart[V+1] = sum(K)+1;
  for (v in 2:V) vstart[v] = vstart[v-1] + K[v-1];
}

parameters {
  real<lower=0> sigma_u[V,D];
  vector[sum(K)] u_raw[D];           //raw random effects
  matrix[P,D] beta;
  vector[D] alpha;
  vector<lower=0>[sum(L)-2*D] cutdiff;
}

transformed parameters {
  vector[sum(L)-D] cutpoint;
  vector[N] eta[D];
  {
    int start = 1;

    for (d in 1:D) {
      
      cutpoint[start]=0;
      for (l in 1:(L[d]-2) ) {
        cutpoint[start+l] = cutpoint[start+l-1] + cutdiff[start+l-d];
      }
      start = start + L[d]-1;
      
      eta[d] = alpha[d] + X*beta[,d]; //+ sigma_g * L_A * g_raw;
      for (v in 1:V) {
        vector[K[v]] u;
        u = sigma_u[v,d] * u_raw[d][vstart[v]:(vstart[v+1]-1)];
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
    
    sigma_u[,d] ~ normal(0,2);
    u_raw[d] ~ normal(0,1);
  }
  to_vector(beta) ~ student_t(4,0,2);
  alpha ~ normal(0,2);
  cutdiff ~ normal(0,5);
}

