n <- 2000
ncat <- rep(3,10)
d <- length(ncat)
K <- 3
pi <- gtools::rdirichlet(1,rep(1,K))
alpha <- rnorm(d)
gamma <- rnorm(d*K)*1.5 %>% matrix(nrow=d,ncol=K)
Z <- LaplacesDemon::rcat(n,as.vector(pi))
cp <- list()
y <- matrix(nrow=n,ncol=d)
for (i in 1:d) {
  cp[[i]] <- c(0,runif(ncat[i]-2)*3) %>% sort()
  y[,i] <- cut(alpha[i] + gamma[i,Z] + rlogis(n),breaks = c(Inf,cp[[i]],-Inf),ordered_result = T)
}

library(rstan)
standat <- list(N=n,D=d,L=ncat,K=K,Y=y)
ordregtop <- stan_model("~/code/OrdRegMix/mordreg_topic.stan")
stanfit <- sampling(ordregtop,data=standat,chains=1,iter=300)
