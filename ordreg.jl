import DataArrays.cut, Lazy.@>, Lazy.@>>, Lazy.@as
using FreqTables.freqtable
using DataFrames, Distributions, LogTopReg


###### prep outer data
n = 100;
nd = repeat([25],inner=[n]);
nobs = sum(nd);
p = 4;
dim = 2;
cp = [[-Inf, 0, 2.0, Inf]];

βout = randn(p,dim);
Xin = randn(p,n);
Xf = repeat(Xin',inner=(nd[1],1));

v = 1;
l = [50]
Xr = [ModelMatrix(ModelFrame(@formula(y~0+ X),DataFrame(y=fill(0,n),
      X=@pdata(repeat(vcat(1:i),inner=div(n,i)))))).m for i=l]
for i = 1:v
    Xr[i] = repeat(Xr[i],inner=(nd[1],1))
end
σ2_u_out = rand(v,dim);
uout = [[sqrt(σ2_u_out[i,d]) .* randn(l[i]) for i=1:v] for d=1:dim];

#### prep inner model
σ_μ = 1.0;
σ = 0.5;
σ_β = 0.0;
K = 3;
αout = rand(dim,K)*4-2;


μ = randn(K)*σ_μ;
βin = randn(K,p)*σ_β;
σy = 1.0;
η = Array{Float64}(K,n);
nk = Array{Int64}(K,n);
for i in 1:n
    η[:,i] = randn(K).*σ .+ μ + βin*Xin[:,i];
    nk[:,i] = rand(Multinomial(nd[i],softmax(η[:,i])));
end
docrng = Vector{UnitRange}(n);
docrng[1] = 1:nd[1];
for i=2:n docrng[i] = (last(docrng[i-1])+1):(nd[i]+last(docrng[i-1])); end
zin = [vcat([fill(k,nk[k,i]) for k=1:K]...) for i=1:n];
zin = vcat(zin...);
#yin = hcat(yin...);

pss0 = VectorPosterior(NormalMeanPosterior(0.0,0.001,σy),2);
hyper = hyperparameter(τ0_u = 0.01,ν0_u=0.33,ν0_σ2η=0.01,σ0_σ2η = 0.01,τ0_τ = 0.01,ν0_τ = 0.01,τ_β = 1e-6);
fit = TLMMfit([init_params(K,length(docrng),size(Xin,1))],
            Matrix{VectorPosterior{NormalMeanPosterior}}(0,0),
            pss0,hyper);
topic = [Normal(randn(),σy) for i=1:dim, j=1:K];

##### back to outer
y = Matrix{Int64}(nobs,dim);
for d in 1:dim
    tmp = cut(αout[d,zin] + Xf*βout[:,d] + Ztu(Xr,uout[d]) + randn(nobs),cp[d]);
    setlevels!(tmp,map(string,1:(length(cp[d])-1)));
    y[:,d] = parse.(tmp);
end

iter = 1000;
prior = Dict{Symbol,Float64}()
prior[:σ2_β] = 5.0;
#prior[:σ2_α] = 5.0;
prior[:τ0_δ] = 0.01;
prior[:ν0_δ] = 0.01;
prior[:τ0_u] = 0.01;
prior[:ν0_u] = 0.01;
γ = [-Inf,0,1,Inf];
l = size.(Xr,1);

samp = Dict{Symbol,AbstractArray}();
samp[:β] = Array{Float64,3}(p,dim,iter);
samp[:α] = Array{Float64,2}(dim,iter);
samp[:δ2] = Array{Float64,2}(dim,iter);
samp[:u] = Array{Vector{Float64},3}(v,dim,iter);
samp[:σ2_u] = Array{Float64,3}(v,dim,iter);

β = randn(p);
δ2 = repeat([1.0],inner=dim);
u = [randn(l[1])];
σ2_u = rand(v)*2;
z = [rand( Truncated( Normal(0,2), cp[d][y[i,d]],cp[d][y[i,d]+1]) ) for i=1:nobs, d=1:dim];

resid = Array{Float64,2}(n,dim);

#initialize topic model
K = 3;
topic = [Normal(randn(),sqrt(δ2[d])) for d=1:dim, j=1:K];

for t in 1:iter

    for d in 1:dim resid[:,d] = z[:,d] - Xf*β[:,d] + Ztu(Xr,u[d])); end
    α, topic, fit = sample_α(resid',Xf,fit,topic,K);

    for d in 1:dim
        resid[:,d] = z[:,d] - (α[:,d] + Ztu(Xr,u[d]));
        β[:,d] = sample_β(resid,Xf,prior[:σ2_β],δ2[d]);

        α_Xβ = α[:,d] + Xf*β[:,d];
        resid[:,d] = z[:,d] - α_Xβ;
        sample_u!(u[d],σ2_u[:,d],resid[:,d],Xr,δ2[d]);
        σ2_u = sample_σ2.(u[d],prior[:ν0_u],prior[:τ0_u],δ2[d]);

        α_Xβ_Zu = α_Xβ + Ztu(Xr,u[d]);
        z[:,d] = sample_z.(α_Xβ_Zu,y[:,d],γ[d],δ2[d]);

        resid[:,d] = z[:,d] - α_Xβ_Zu;
        δ2[d] = sample_σ2(resid[:,d],prior[:ν0_δ],prior[:τ0_δ]);

        for k in 1:K, topic[d,k] = Normal(topic[d,k].μ,sqrt(δ2[d])); end

    end

    samp[:β][:,:,t] = β;
    samp[:α][:,t] = α;
    samp[:δ2][:,t] = δ2;
    samp[:u][:,:,t] = deepcopy(u);
    samp[:σ2_u][:,:,t] = σ2_u;
end


function Ztu(Z,u)
    return reduce(+,Z.*u)
end
