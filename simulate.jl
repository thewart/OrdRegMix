using DataFrames, Distributions, LogTopReg, CategoricalArrays,
    StatsModels, StatsBase, StatsFuns
import LogTopReg.init_params, LogTopReg.init_params!, StatsBase.cov2cor
include("/home/seth/code/OrdRegMix/hybridout.jl");
include("/home/seth/code/OrdRegMix/samplers.jl");

cov2cor(C) = cov2cor(C,sqrt.(diag(C)));

###### prep outer data
n = 400;
nd = repeat([20],inner=[n]);
docrng = nd_to_docrng(nd);
nobs = sum(nd);
p = 3;
cp = [-Inf, 0, 1.0, Inf];
dim = 2;

βout = randn(p,dim)*0.5;
#βout = zeros(p,dim);
Xin = randn(n,p);
#Xin = vcat(ones(1,n),randn(p-1,n));
Xf = vcat(map(i-> repeat(Xin[i:i,:],outer=(nd[i],1)),1:n)...);

v = 2;
l = [200,100];
Xr = [ModelMatrix(ModelFrame(@formula(y~0+ X),DataFrame(y=fill(0,nobs),
      X=CategoricalArray(shuffle(repeat(collect(1:i),inner=div(nobs,i))))))).m
      for i in l];
Xr = [convert.(Bool,i) for i in Xr];
#Xr = randn(size(Xr));

#σ2_u_out = rand(v,dim);
σ2_u_out = [rand(InverseWishart(dim*2,eye(dim))) for i in 1:v];
#σ2_u_out = cov2cor.(σ2_u_out).*0.5;
#σ2_u_out = [diagm(rand(dim)) for i in 1:v];
uout = [ (chol(σ2_u_out[i])*randn(dim,l[i]))' for i in 1:v];
#σ2_out = rand(dim);
σ2_out = fill(1.0,dim);
#### prep inner model
#σ_μ = 1.0;
σ = 0.5;
σ_β = 0.0;
K = 2;
#αout = randn(dim,K)*0.75;
αout = [1 -1;-1 1].*1.5;
#αout = rand(dim,K)*0;

#μ = randn(K)*σ_μ;
μ = fill(0.0,K);
βin = randn(K,p)*σ_β;
η = Array{Float64}(K,n);
nk = Array{Int64}(K,n);
for i in 1:n
    η[:,i] = randn(K).*σ .+ μ + βin*Xin[i,:];
    nk[:,i] = rand(Multinomial(nd[i],softmax(η[:,i])));
end
docrng = Vector{UnitRange}(n);
docrng[1] = 1:nd[1];
for i=2:n docrng[i] = (last(docrng[i-1])+1):(nd[i]+last(docrng[i-1])); end
zin = [vcat([fill(k,nk[k,i]) for k=1:K]...) for i=1:n];
zin = vcat(zin...);
#yin = hcat(yin...);

# pss0 = [NormalMeanPosterior(0.0,0.001,1.0) for i=1:dim]; #note: 1.0 value is a placeholder
# hyper = hyperparameter(τ0_u = 0.01,ν0_u=0.33,ν0_σ2η=0.01,σ0_σ2η = 0.01,τ0_τ = 0.01,ν0_τ = 0.01,τ_β = 1e-6);
# fit = TLMMfit([init_params(K,length(docrng),size(Xin,1))],
#             Matrix{Vector{NormalMeanPosterior}}(0,0),
#             pss0,hyperparameter());

##### back to outer
y = Matrix{Int64}(nobs,dim);
eta = αout[:,zin]' + Xf*βout + Ztu(Xr,uout) + randn(nobs,dim).*sqrt.(σ2_out)';
for d in 1:dim
    tmp = cut(eta[:,d],cp,labels=map(string,1:(length(cp)-1)));
    y[:,d] = parse.(tmp);
    #y[:,d] = αout[d,zin] + Xf*βout[:,d] + Ztu(Xr,uout[d]) + randn(nobs)*sqrt(σ2_out[d]);
end

#### initialize
#γ = [-Inf,0,1,Inf];
#l = size.(Xr,1);
