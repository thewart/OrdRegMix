using DataFrames, Distributions, LogTopReg, CategoricalArrays,
    StatsModels, StatsBase, StatsFuns

###### prep outer data
n = 200;
nd = repeat([20],inner=[n]);
docrng = nd_to_docrng(nd);
nobs = sum(nd);
p = 3;
cp = [-Inf, 0, 1.0, Inf];
dim = 4;

#βout = randn(p,dim);
βout = zeros(p,dim);
Xin = randn(n,p);
#Xin = vcat(ones(1,n),randn(p-1,n));
Xf = vcat(map(i-> repeat(Xin[i:i,:],outer=(nd[i],1)),1:n)...);

v = 1;
l = 100;
Xr = ModelMatrix(ModelFrame(@formula(y~0+ X),DataFrame(y=fill(0,n),
      X=CategoricalArray(repeat(vcat(1:l),inner=div(n,l)))))).m;
Xr = [convert.(Bool,repeat(Xr,inner=(nd[1],1)))];
#Xr = randn(size(Xr));

σ2_u_out = rand(v,dim);
#σ2_u_out = zeros(v,dim);
uout = [hcat([sqrt(σ2_u_out[d]) .* randn(l) for d=1:dim]...)];
σ2_out = rand(dim);

#### prep inner model
σ_μ = 1.0;
σ = 0.5;
σ_β = 0.0;
K = 2;
αout = randn(dim,K)*0.75;
#αout = rand(dim,K)*0;

μ = randn(K)*σ_μ;
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
eta = αout[:,zin]' + Xf*βout + Ztu(Xr,uout) + randn(nobs,dim).*sqrt.(σ2_out)'
for d in 1:dim
    tmp = cut(eta[:,d],cp,labels=map(string,1:(length(cp)-1)));
    y[:,d] = parse.(tmp);
    #y[:,d] = αout[d,zin] + Xf*βout[:,d] + Ztu(Xr,uout[d]) + randn(nobs)*sqrt(σ2_out[d]);
end

#### initialize
#γ = [-Inf,0,1,Inf];
#l = size.(Xr,1);
