library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

ordreg <- stan_model("~/code/OrdRegMix/ordreg.stan")

n <- 1000
p <- 3
ncat <- 5
cp <- c(0,runif(ncat-2)*3) %>% sort()
beta <- rnorm(p)
alpha <- rnorm(1)
X <- rnorm(n*p) %>% matrix(nrow=n)
y <- cut(alpha + X %*% beta + rlogis(n),breaks = c(Inf,cp,-Inf),labels=F)

standat <- list(y=y,X=X,N=n,P=p,L=ncat)

stanfit <- sampling(ordreg,standat,chains=1,iter=500,warmup=100,init=0)


mordreg <- stan_model("~/code/OrdRegMix/mordreg.stan")

n <- 1000
p <- 3
ncat <- c(2,3,4)
d <- length(ncat)

beta <- rnorm(p*d) %>% matrix(nrow=p)
alpha <- rnorm(d)
X <- rnorm(n*p) %>% matrix(nrow=n)
y <- matrix(nrow=n,ncol=d)
cp <- list()
for (i in 1:d) {
  cp[[i]] <- c(0,runif(ncat[i]-2)*3) %>% sort()
  y[,i] <- cut(alpha[i] + X %*% beta[,i] + rlogis(n),breaks = c(Inf,cp[[i]],-Inf),labels=F)
}
standat <- list(Y=y,X=X,N=n,P=p,D=d,L=ncat)
stanfit <- sampling(mordreg,standat,chains=1,iter=500,warmup=250,init=0)

standat <- list(y=y[,1],X=X,N=n,P=p,L=ncat[1])
stanfit <- sampling(ordreg,standat,chains=1,iter=500,warmup=250,init=0)
