addprocs(3);
path = "/home/seth/analysis/OrdRegMix/120518fF_rank/";
include("/home/seth/code/OrdRegMix/readindata.jl")

#read in initialization
Kvec = [1,2,4,6,8];
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

    # iter = 10000;
    # thin = 20;

    hy = hyu = hyperparameter(ν0_τ=1.0,τ0_τ=1.0,τ_β=5.0);
    hy[:a] = 1.0;
    hy[:ν0] = 2.0;
    hyin = hyperparameter(σ0_σ2η=0.1,ν0_σ2η=1.0);

    hyinu = hyin0 = hyperparameter(σ0_σ2η=1e-12,ν0_σ2η=1e12);

    hy0 = hyperparameter(ν0_τ=1.0,τ0_τ=1.0,τ_β=5.0);
    hy0[:a] = 1e-12;
    hy0[:ν0] = 1e12;

    hyvec = [hy,hyu,hy0];
    hyinvec = [hyin,hyinu,hyin0];

    f1 = (hy,hyin) -> lmmtopic(Y,Xf,Xr,Xin,docrng,K,init=init,hy=hy,hyin=hyin,
        initin=initin,iter=10000,thin=20);
    f2 = fit -> lppd(Y,Xf,Xr,docrng,fit);

    @time f1out = pmap(f1,hyvec,hyinvec);
    f2out = pmap(f2,f1out);
    foof[i] = f1out[1];
    foofu[i] = f1out[2];
    foof0[i] = f1out[3];
    lp[i] = f2out[1];
    lpu[i] = f2out[2];
    lp0[i] = f2out[3];
end

@save string(path,"topicfits.jld") foof foofu foof0 lp lpu lp0

@load string(path,"topicfits.jld") foof;

d = size(Y,2);
n = size(Xr[1],2);
i = 5;
fit = foof[i];
py = Array{Float64}(2^d,n,length(fit)));
for j in eachindex(fit)
    s = fit[j];
    for l in 1:n
        pi = softmax(s.tlmm.η[:,l]);
        u = s.u[1][l,:];
        pyd = calc_pyd(s.α,u);
        py[:,l,j] = calc_py(pyd,pi);
    end
end
