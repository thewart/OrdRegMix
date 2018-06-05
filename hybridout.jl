mutable struct HYBRIDsample
    α::Matrix{Float64}
    β::Matrix{Float64}
    u::Vector{Matrix{Float64}}

    σ2::Vector{Float64}
    Σ::Vector{Matrix{Float64}}

    tlmm::TLMMsample

    HYBRIDsample() = new();
end

function init_params!(samp::HYBRIDsample,p::Int,l::Vector{Int},dim::Int,K::Int)
    v = length(l);
    samp.β = randn(p,dim);
    samp.σ2 = rand(dim);
    samp.u = [randn(i,dim) for i=l];
    samp.Σ = [(samp.u[i]'samp.u[i])./l[i] for i=1:v];
    samp.α = randn(dim,K);
    #samp.u = [randn(i) for i=l, d=1:dim];
    #samp.σ2_u = rand(v,dim)*2;

    return samp
end

function init_params(p::Int,l::Vector{Int},dim::Int,K::Int)
    return init_params!(HYBRIDsample(),p,l,dim,K)
end

function nd_to_docrng(nd)
    bounds = vcat(0,cumsum(nd));
    n = length(nd);
    return [(bounds[i]+1):bounds[i+1] for i=1:n];
end

function lmmtopic(y,Xf,Xr,Xin,docrng,K;hy=hyperparameter(),hyin=hyperparameter(),
    init=init_params(size(Xf,2),size.(Xr,2),size(y,2),K),
    initin=init_params(K,length(docrng),size(Xin,2)),
    cp=fill(1,size(y,2)),iter=1000,thin=1)

    p = size(Xf,2);
    n,dim = size(y);
    v = length(Xr);
    l = size.(Xr,2);

    ν0_α = 0.001;
    μ0_α = 0.0;

    if !isdefined(initin,:z)
        initin.z = init_z(initin.η,K,docrng);
    end

    fit = TLMMfit([initin],
        Matrix{Vector{NormalMeanPosterior}}(0,0),
        [NormalMeanPosterior(μ0_α,ν0_α,1.0) for i=1:dim],
        hyin);

    #s = init_params(p,l,dim);
    s = deepcopy(init);

    resid = Array{Float64,2}(n,dim);
    α_expand =  s.α[:,initin.z]';
    η = α_expand + Xf*s.β + Ztu(Xr,s.u);
    z = [sample_z(η[i,d],y[i,d],s.σ2[d],cp[d]) for i=1:n, d=1:dim];

    saveiter = thin:thin:iter;
    nsave = length(saveiter);
    iter = maximum(saveiter);

    samples = Vector{HYBRIDsample}(nsave);

    iLu = Vector{Matrix{Float64}}(v);

    ##### run model
    for t in 1:iter

        resid = z - Xf*s.β - Ztu(Xr,s.u);
        for d in 1:dim
            fit.prior[d].σ2 = s.σ2[d];
        end
        s.α, fit = sample_α(resid',Xin',fit,docrng,K);
        α_expand =  s.α[:,fit.θ[1].z]';
        resid = z - α_expand + Xf*s.β;

        if v==1
            sample_u!(s.u[1],resid,Xr[1],
                diagm(sqrt.(s.σ2))*s.Σ[1]*diagm(sqrt.(s.σ2)),s.σ2);
            s.Σ[1] = sample_Σ(s.u[1].*inv.(sqrt.(s.σ2))',
                diag(s.Σ[1]),hy[:ν],hy[:a0]);
            iLu[1] = s.u[1]*chol(Hermitian(inv(s.Σ[1])))';
        else
            for i in 1:v
                noti = setdiff(1:v,i);
                sample_u!(s.u[i],resid-Ztu(Xr[noti],s.u[noti]),Xr[i],
                    diagm(sqrt.(s.σ2))*s.Σ[i]*diagm(sqrt.(s.σ2)),s.σ2);
                s.Σ[i] = sample_Σ(s.u[i].*inv.(sqrt.(s.σ2))',
                    diag(s.Σ[i]),hy[:ν],hy[:a0]);
                iLu[i] = s.u[i]*chol(Hermitian(inv(s.Σ[i])))';
            end
        end
        α_Zu = α_expand + Ztu(Xr,s.u);
        resid = z - α_Zu;

        for d in 1:dim
            s.β[:,d] = sample_β(resid[:,d],Xf,hy[:τ_β]*s.σ2[d],s.σ2[d]);

            α_Xβ_Zu = α_Zu[:,d] + Xf*s.β[:,d];
            z[:,d] = sample_z.(α_Xβ_Zu,y[:,d],s.σ2[d],cp[d]);

            resid[:,d] = z[:,d] - α_Xβ_Zu;
            bigvec = vcat(resid[:,d],(s.α[d,:].-μ0_α)*sqrt(ν0_α),
                    s.β[:,d]./sqrt(hy[:τ_β]),
                    vcat([iLu[i][:,d] for i in 1:v]...));
            s.σ2[d] = sample_σ2(bigvec,hy[:ν0_τ],hy[:τ0_τ]);
        end

        if t ∈ saveiter
            j = findin(saveiter,t)[1];
            samples[j] = deepcopy(s);
            samples[j].tlmm = deepcopy(fit.θ[1]);
        end
    end

    return samples
end

import LogTopReg.gf
function gf(value::Vector{HYBRIDsample},name::Symbol)
    nd = ndims(getfield(value[1],name));
    return cat(nd+1,getfield.(value,name)...)
end

import Base.mean
function mean(samps::Vector{HYBRIDsample})
    meanfit = HYBRIDsample();

    tmp = gf(samps,:α);
    meanfit.α = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));
    tmp = gf(samps,:β);
    meanfit.β = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));
    tmp = gf(samps,:σ2);
    meanfit.σ2 = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));
    tmp = gf(samps,:σ2_u);
    meanfit.σ2_u = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));

    v,dim = size(meanfit.σ2_u,1,2)
    tmp = gf(samps,:u);
    meanfit.u = [squeeze(mean(hcat(tmp[i,j,:]...),2),2) for i in 1:v, j in 1:dim];

    meanfit.tlmm = mean([i.tlmm for i in samps]);

    return meanfit
end

###### prediction for new samples from old groups
# function probord_pdf(y::Int64,η::Float64,cp)
#     dist = Normal(η,1);
#     return log(cdf(dist,cp[y+1]) - cdf(dist,cp[y]));
# end
#
# function probordmix_pdf(y,η,α,π,cp);
#     K = length(π);
#     dim = length(η);
#     lp = Vector{Float64}(length(π));
#     for k ∈ 1:K
#         lp[k] = log(π[k]);
#         for d ∈ 1:dim
#             lp[k] += probord_pdf(y[d],η[d] + α[d,k],cp[:,d])
#         end
#     end
#     return logsumexp(lp)
# end
#
# import LogTopReg.lppd
# function lppd(y,Xf,Xr,docrng,cp,s::HYBRIDsample)
#     n = length(docrng);
#     lp = zeros(n);
#     dim = length(cp);
#     cpf = repeat([-Inf,0,1,Inf],outer=(1,dim));
#     cpf[3,:] = cp;
#
#     η = Matrix{Float64}(size(y));
#     for d ∈ 1:dim
#         η[:,d] = Xf*s.β[:,d] + Ztu(Xr,s.u[:,d]);
#     end
#
#     for i ∈ 1:n
#         irng = docrng[i];
#         if size(s.tlmm.η,2) == 1
#             π = softmax(s.tlmm.η[:,1]);
#         else
#             π = softmax(s.tlmm.η[:,i]);
#         end
#
#         for j ∈ irng
#             lp[i] += probordmix_pdf(y[j,:],η[j,:],s.α,π,cpf);
#         end
#     end
#
#     return lp;
# end
#
# function lppd(y,Xf,Xr,docrng,cp,fit::Vector{HYBRIDsample})
#     iter = length(fit);
#     n = length(docrng);
#     lp = Matrix{Float64}(n,iter);
#
#     for i ∈ 1:iter
#         lp[:,i] = lppd(y,Xf,Xr,docrng,cp,fit[i]);
#     end
#
#     return lp
# end
#
# function splitfolds(n,k)
#     nrep = ceil(Int,n/k);
#     return shuffle!(repeat(1:k,inner=nrep))[1:n]
# end
#
# function foldup(docrng,k)
#     folds = splitfolds.(length.(docrng),k);
#     flatfolds = vcat(folds...);
#     foldrng = mapslices(nd_to_docrng,[countnz(juh[i].==j) for i=1:length(docrng), j=1:k],1);
#     return flatfolds, foldrng
# end
