library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

ordreg <- stan_model("~/code/OrdRegMix/ordreg.stan")

r <- 500
rep <- 1
p <- 2
ncat <- 3
cp <- c(0,runif(ncat-2)*3) %>% sort()
K <- c(5,10)
sigma_u <- c(.75,0.1)
#sigma_g <- 0.0
beta <- rnorm(p)
alpha <- rnorm(1)
X <- rnorm(r*p) %>% matrix(nrow=r)
#L_A <- chol(MCMCpack::riwish(n+1,diag(r))) %>% t()
u <- lapply(1:length(K),function(x) sigma_u[x] * rnorm(K[x]))
Z <- sapply(K,function(x) sample(1:x,r,replace=T))
eta <- alpha + X %*% beta #+ sigma_g * L_A %*% rnorm(r)
for (k in 1:length(K)) eta <- eta + u[[k]][Z[,k]]

y <- matrix(0,nrow=r,ncol=ncat)
for (i in 1:r) y[i,] <- cut(eta[i] + rlogis(10),breaks = c(Inf,cp,-Inf),ordered_result = T) %>% table()

standat <- list(y=y,X=X,N=r,P=p,V=length(K),L=ncat,K=K,Z=Z)

stanfit <- sampling(ordreg,standat,chains=1,iter=500,warmup=200,pars=c("alpha","beta","sigma_u","cutpoint"))



mordreg <- stan_model("~/code/OrdRegMix/mordreg.stan")

n <- 500
p <- 3
ncat <- c(2,3)
d <- length(ncat)

beta <- rnorm(p*d) %>% matrix(nrow=p)
alpha <- rnorm(d)
X <- rnorm(n*p) %>% matrix(nrow=n)
K <- c(5,10)
Z <- sapply(K,function(x) sample(1:x,n,replace=T))
sigma_u <- matrix(c(1.0,0.1,0.5,1.5),nrow = 2)

cp <- list()
ytmp <- list()
for (i in 1:d) {
  cp[[i]] <- c(0,runif(ncat[i]-2)*3) %>% sort()
  u <- lapply(1:length(K),function(x) sigma_u[x,i] * rnorm(K[x]))
  eta <- alpha[i] + X %*% beta[,i]
  for (k in 1:length(K)) eta <- eta + u[[k]][Z[,k]]
  
  ytmp[[i]] <- matrix(nrow=n,ncol=ncat[i])
  for (j in 1:n) ytmp[[i]][j,] <- cut(eta[j] + rlogis(10),breaks = c(Inf,cp[[i]],-Inf),ordered_result = T) %>% table()
}
y <- do.call(cbind,args=ytmp)
standat <- list(Y=y,X=X,N=n,P=p,D=d,L=ncat,V=length(K),K=K,Z=Z)
stanfit <- sampling(mordreg,standat,chains=2,iter=500,warmup=250,pars=c("alpha","beta","sigma_u","cutpoint"))


standat <- list(y=y[,3:5],X=X,N=n,P=p,V=length(K),L=ncat[2],K=K,Z=Z)
stanfit <- sampling(ordreg,standat,chains=1,iter=400,warmup=200,pars=c("alpha","beta","sigma_u","cutpoint"))


Z <- matrix(1,n,length(K))
ynew <- matrix(0,n,4)
ynew[cbind(1:n,y[,3])] <- 1
standat <- list(y=ynew,X=X,N=n,P=p,V=length(K),L=ncat[3],K=K,Z=Z)
stanfit <- sampling(ordreg,standat,chains=1,iter=400,warmup=200,pars=c("alpha","beta","sigma_u","cutpoint"))
