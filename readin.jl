@everywhere using DataFrames, Distributions, LogTopReg, StatsBase, StatsFuns, Distributions
import LogTopReg.init_params, LogTopReg.init_params!
@everywhere include("/home/seth/code/OrdRegMix/hybridout.jl");
include("/home/seth/code/OrdRegMix/samplers.jl");
@everywhere include("/home/seth/code/OrdRegMix/entropy.jl")

#read in data
path = "/home/seth/analysis/OrdRegMix/071018fF/";
v = sum(getindex.(readdir(path),[1:2]) .== "Xr");
Xf = readcsv(string(path,"Xf.csv"));
Xr = [readcsv(string(path,"Xr",i,".csv")) for i=1:v];
Xr = [convert.(Bool,X) for X in Xr];
Y = readcsv(string(path,"Y.csv"),Int64);
nd = readcsv(string(path,"docrng.csv"),Int64)[:];
docrng = nd_to_docrng(nd);
Xin = zeros(length(docrng),1);
# Xin = readcsv(string(path,"Xin.csv"));

fv = false;

# nfolds = 4;
# folds = foldup(docrng,nfolds);

#read in initialization
Kvec = collect(2:8);
foof = Vector{Vector{HYBRIDsample}}(length(Kvec));
lp = Vector{Matrix{Float64}}(length(Kvec));
# cv = Matrix{Matrix{Float64}}(nfolds,length(Kvec));

for i in 1:length(Kvec)
    K = Kvec[i];
    suff = string("_k",K,".csv")
    init = HYBRIDsample();
    #init = init_params(size(Xf,2),size.(Xr,2),size(Y,2));
    init.σ2 = vec(readcsv(string(path,"sigma_0",suff)).^2);
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

    resp = readcsv(string(path,"r_0",suff));
    initin = init_params(size(init.α,2),length(docrng),size(Xin,2));
    iz = vec(rand.(mapslices(Categorical,resp,2)));
    initin.z = iz;

    if fv
        init.σ2[:] = 1.0;
        cp = inv.(vec(readcsv(string(path,"sigma_0",suff))));
    else
        cp = ones(size(Y,2));
    end

    hy = hyperparameter(ν0_τ=1.0,τ0_τ=1.0,τ_β=5.0);
    hy[:a] = 1.0;
    hy[:ν0] = 2.0;
    #hack observers
    init.σ2_u[2] = diagm(diag(init.σ2_u[2]));

    # @time foof[i] = lmmtopic(Y,Xf,Xr,Xin,docrng,init=init,hy=hy,
    # initin=initin,K,iter=10,thin=1,cp=cp,fixedvar=fv);
    # @time lp[i] = lppd(Y,Xf,Xr,docrng,cp,foof[i]);
    # fit = mean(foof[i]);
    # @time H[i] = calc_cond_entropy_parallel_2(fit,cp,2);

    @time foof[i] = lmmtopic(Y,Xf,Xr,Xin,docrng,init=init,hy=hy,
    initin=initin,K,iter=20000,thin=40,cp=cp,fixedvar=fv);
    @time lp[i] = lppd(Y,Xf,Xr,docrng,cp,foof[i]);

    # for j in 1:nfolds
    #     trainfolds = setdiff(1:nfolds,j);
    #     train = [[g ∈ trainfolds for g in f] for f in folds];
    #     test = [f .== j for f in folds];
    #     trainrng =  nd_to_docrng([countnz(m) for m in train]);
    #     testrng = nd_to_docrng([countnz(m) for m in test]);
    #
    #     train = vcat(train...);
    #     test = vcat(test...);
    #
    #     Xr_train = [X[train,:] for X in Xr];
    #     Xr_test = [X[test,:] for X in Xr];
    #     initin.z = iz[train];
    #
    #     @time fit = lmmtopic(Y[train,:],Xf[train,:],Xr_train,Xin,trainrng,init=init,hy=hy,
    #         initin=initin,K,iter=5000,thin=5,cp=cp,fixedvar=fv);
    #     cv[j,i] = lppd(Y[test,:],Xf[test,:],Xr_test,testrng,cp,fit);
    # end

end

rng = 201:3:500;
H = Matrix{Tuple{Float64,Float64,Float64}}(length(rng),length(Kvec));
for i in 1:length(Kvec)
    @time H[:,i] = pmap(calc_cond_entropy_2,foof[i][rng],fill(150,length(rng)));
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
