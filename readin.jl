using DataFrames, Distributions, StatsBase, StatsFuns
include("/home/seth/code/LogisticTopicRegression/LogTopReg.jl");
using LogTopReg
import LogTopReg.init_params, .LogTopReg.init_params!
include("/home/seth/code/OrdRegMix/hybridout.jl");
include("/home/seth/code/OrdRegMix/samplers.jl");
include("/home/seth/code/OrdRegMix/entropy.jl");

#read in data
path = "/home/seth/analysis/OrdRegMix/111418fF_kin/";
v = sum(getindex.(readdir(path),[1:2]) .== "Xr");
Xf = readcsv(string(path,"Xf.csv"));
Xr = [readcsv(string(path,"Xr",i,".csv")) for i=1:v];
Xr = [convert.(Bool,X) for X in Xr];
Y = Bool.(readcsv(string(path,"Y.csv"),Int64));
nd = readcsv(string(path,"docrng.csv"),Int64)[:];
docrng = nd_to_docrng(nd);
Xin = zeros(length(docrng),1);
# Xin = readcsv(string(path,"Xin.csv"));

# nfolds = 4;
# folds = foldup(docrng,nfolds);

#read in initialization
Kvec = collect(1:2:6);
foof = Vector{Vector{HYBRIDsample}}(length(Kvec));
lp = Vector{Matrix{Float64}}(length(Kvec));
foofu = Vector{Vector{HYBRIDsample}}(length(Kvec));
lpu = Vector{Matrix{Float64}}(length(Kvec));
foof0 = Vector{Vector{HYBRIDsample}}(length(Kvec));
lp0 = Vector{Matrix{Float64}}(length(Kvec));

# cv = Matrix{Matrix{Float64}}(nfolds,length(Kvec));

for i in 1:length(Kvec)
    K = Kvec[i];
    suff = string("_k",K,".csv")
    init = HYBRIDsample();
    #init = init_params(size(Xf,2),size.(Xr,2),size(Y,2));
    init.u = [readcsv(string(path,"ranef",i,"_0",suff)) for i=1:v];
    #init.u = [uinit[i][:,d] for i=1:v,d=1:size(Y,2)];
    init.α = readcsv(string(path,"alpha_0",suff));
    init.β = readcsv(string(path,"fixef_0",suff));
    #init.σ2_u = var.(init.u./sqrt.(init.σ2'));
    init.σ2_u = cov.(init.u);
    for Σ in init.σ2_u
        for j in 1:size(Σ,1)
            if Σ[j,j]<0.01 Σ[j,j] = 0.01; end
        end
    end

    #init.σ2_u[init.σ2_u.<0.01] = 0.01;
    #K = size(init.α,2);
    if K==1
        iz = ones(size(Y,1));
    else
        resp = readcsv(string(path,"r_0",suff));
        iz = vec(rand.(mapslices(Categorical,resp,2)));
    end
    initin = init_params(size(init.α,2),length(docrng),size(Xin,2));
    initin.z = iz;
    initin.η = exp.([mean(initin.z[d].==k) for d in docrng, k in 1:K])';
    initin.μ = mean(initin.η,2)[:,1];

    hy = hyperparameter(ν0_τ=1.0,τ0_τ=1.0,τ_β=5.0);
    hy[:a] = 1.0;
    hy[:ν0] = 2.0;

    iter = 10000;
    thin = 10;

    hyin = hyperparameter(σ0_σ2η=0.1,ν0_σ2η=1.0);
    @time foof[i] = lmmtopic(Y,Xf,Xr,Xin,docrng,K,init=init,hy=hy,
        initin=initin,iter=iter,thin=thin);
    @time lp[i] = lppd(Y,Xf,Xr,docrng,foof[i]);

    hyin = hyperparameter(σ0_σ2η=1e-12,ν0_σ2η=1e12);
    @time foofu[i] = lmmtopic(Y,Xf,Xr,Xin,docrng,K,init=init,hy=hy,hyin=hyin,
        initin=initin,iter=iter,thin=thin);
    @time lpu[i] = lppd(Y,Xf,Xr,docrng,foofu[i]);

    hy[:a] = 1e-12;
    hy[:ν0] = 1e12;
    @time foof0[i] = lmmtopic(Y,Xf,Xr,Xin,docrng,K,init=init,hy=hy,hyin=hyin,
        initin=initin,iter=iter,thin=thin);
    @time lp0[i] = lppd(Y,Xf,Xr,docrng,foof0[i]);
end

rng = 1:2:500;
H = Matrix{Tuple{Float64,Float64,Float64}}(length(rng),length(Kvec)-1);
H0 = Matrix{Tuple{Float64,Float64,Float64}}(length(rng),length(Kvec)-1);
for i in 1:(length(Kvec)-1)
    # @time H0[:,i] = pmap(calc_cond_entropy,foof[i+1][rng],fill(100,length(rng)),fill(true,length(rng)));
    @time H[:,i] = pmap(calc_cond_entropy,foof[i+1][rng],fill(100,length(rng)));
end

initin = init_params(size(init.α,2),1,1);
initin.z = vec(rand.(mapslices(Categorical,resp,2)));
@time foof0 = lmmtopic_fixedvar(Y,Xf,Xr,zeros(1,1),[1:size(Y,1)],init=init,
        initin=initin,K,iter=5000,thin=5,cp=cp);
lp0 = lppd(Y,Xf,Xr,docrng,cp,foof0);

pii = mapslices.(softmax,[foof[i].tlmm.η for i=1:iter],1);
mean(cat(3,pii[501:iter]...),3)

function cH(H)
    return (H[1]-H[2])/H[1]
end

function cHpi(H)
    return (H[1]-H[3])/(H[1]-H[2])
end


function suminf(x)
    bad = isinf.(x);
    return sum(x[!bad])
end
