mutable struct HYBRIDsample
    α::Matrix{Float64}
    β::Matrix{Float64}
    u::Vector{Matrix{Float64}}

    σ2::Vector{Float64}
    σ2_u::Vector{Matrix{Float64}}

    tlmm::TLMMsample

    HYBRIDsample() = new();
end

function init_params!(samp::HYBRIDsample,p::Int,l::Vector{Int},dim::Int)
    v = length(l);
    samp.β = randn(p,dim);
    samp.σ2 = rand(dim);
    samp.u = [randn(i,dim) for i in l];
    samp.σ2_u = [samp.u[i]'samp.u[i]./l[i] for i in 1:v];

    return samp
end

function init_params(p::Int,l::Vector{Int},dim::Int)
    return init_params!(HYBRIDsample(),p,l,dim)
end

function nd_to_docrng(nd)
    bounds = vcat(0,cumsum(nd));
    n = length(nd);
    return [(bounds[i]+1):bounds[i+1] for i=1:n];
end

function lmmtopic(y,Xf,Xr,Xin,docrng,K;hy=hyperparameter(),hyin=hyperparameter(),
    init=init_params(size(Xf,2),size.(Xr,2),size(y,2)),
    initin=init_params(K,length(docrng),size(Xin,2)),
    cp=fill(1,size(y,2)),iter=1000,thin=1,fixedvar=false)

    p = size(Xf,2);
    n,dim = size(y);
    v = length(Xr);
    l = size.(Xr,2);

    ν0_α = 0.1;
    μ0_α = 0.0;
    #s = init_params(p,l,dim);
    s = deepcopy(init);

    if isdefined(initin,:z)
        α_expand =  s.α[:,initin.z]';
        η = Matrix{Float64}(size(y));
        η = α_expand + Xf*s.β + Ztu(Xr,s.u);
        z = [sample_z(η[i,d],y[i,d],s.σ2[d],cp[d]) for i=1:n, d=1:dim];
    else
        z = [sample_z(0.0,y[i,d],s.σ2[d],cp[d]) for i=1:n, d=1:dim];
    end

    fit = TLMMfit([initin],
        Matrix{Vector{NormalMeanPosterior}}(0,0),
        [NormalMeanPosterior(μ0_α,ν0_α,1.0) for i=1:dim],
        hyin);
    resid = Array{Float64,2}(n,dim);

    saveiter = thin:thin:iter;
    nsave = length(saveiter);
    iter = maximum(saveiter);

    samples = Vector{HYBRIDsample}(nsave);

    ##### run model
    for t in 1:iter

        resid = z - Xf*s.β - Ztu(Xr,s.u);
        for d in 1:dim
            fit.prior[d].σ2 = s.σ2[d];
        end
        s.α, fit = sample_α(resid',Xin',fit,docrng,K);
        α_expand =  s.α[:,fit.θ[1].z]';

        resid = z - (α_expand + Ztu(Xr,s.u));
        for d in 1:dim
            s.β[:,d] = sample_β(resid[:,d],Xf,hy[:τ_β]*s.σ2[d],s.σ2[d]);
        end

        α_Xβ = α_expand + Xf*s.β;
        resid = z - α_Xβ;
        sample_u!(s.u,s.σ2_u,resid,Xr,s.σ2);

        for j in 1:v
            s.σ2_u[j] = sample_Σ(s.u[j],#*Diagonal(inv.(sqrt.(s.σ2))),
                        diag(s.σ2_u[j]),hy[:ν0],hy[:a]);
        end

        α_Xβ_Zu = α_Xβ + Ztu(Xr,s.u);

        for d in 1:dim
            z[:,d] .= sample_z.(α_Xβ_Zu[:,d],y[:,d],s.σ2[d],cp[d]);
        end

        if !fixedvar
            resid = z - α_Xβ_Zu;
            bigmat = vcat(resid,
                    (s.α'-μ0_α)*sqrt(ν0_α),
                    s.β./sqrt(hy[:τ_β]));#,
                    # vcat([s.u[i]*inv(chol(s.σ2_u[i]))' for i in 1:v]...));
            for d in 1:dim
                s.σ2[d] = sample_σ2(bigmat[:,d],hy[:ν0_τ],hy[:τ0_τ]);
            end
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
    tmp = [cat(3,[i.σ2_u[j] for i in samps]...) for j in 1:2];
    meanfit.σ2_u = [squeeze(mean(i,ndims(i)),ndims(i)) for i in tmp];
    tmp = [cat(3,[i.u[j] for i in samps]...) for j in 1:2];
    meanfit.u = [squeeze(mean(i,ndims(i)),ndims(i)) for i in tmp];

    meanfit.tlmm = mean([i.tlmm for i in samps]);

    return meanfit
end

##### prediction for new samples from old groups
function probord_pdf(y::Int64,η::Float64,σ2,cp)
    dist = Normal(η,sqrt(σ2));
    return log(cdf(dist,cp[y+1]) - cdf(dist,cp[y]));
end

function probordmix_pdf(y,η,α,σ2,π,cp);
    K = length(π);
    dim = length(η);
    lp = Vector{Float64}(length(π));
    for k ∈ 1:K
        lp[k] = log(π[k]);
        for d ∈ 1:dim
            lp[k] += probord_pdf(y[d],η[d] + α[d,k],σ2[d],cp[:,d])
        end
    end
    return logsumexp(lp)
end

import LogTopReg.lppd
function lppd(y,Xf,Xr,docrng,cp,s::HYBRIDsample)
    n = length(docrng);
    lp = zeros(size(y,1));
    dim = length(cp);
    cpf = repeat([-Inf,0,1,Inf],outer=(1,dim));
    cpf[3,:] = cp;

    η = Xf*s.β + Ztu(Xr,s.u);

    for i ∈ 1:n
        irng = docrng[i];
        if size(s.tlmm.η,2) == 1
            π = softmax(s.tlmm.η[:,1]);
        else
            π = softmax(s.tlmm.η[:,i]);
        end

        for j ∈ irng
            lp[j] += probordmix_pdf(y[j,:],η[j,:],s.α,s.σ2,π,cpf);
        end
    end

    return lp;
end

function lppd(y,Xf,Xr,docrng,cp,fit::Vector{HYBRIDsample})
    iter = length(fit);
    # n = size(docrng);
    lp = Matrix{Float64}(size(y,1),iter);

    for i ∈ 1:iter
        lp[:,i] = lppd(y,Xf,Xr,docrng,cp,fit[i]);
    end

    return lp
end

function splitfolds(n,k)
    nrep = ceil(Int,n/k);
    return shuffle!(repeat(1:k,inner=nrep))[1:n]
end

function foldup(docrng,k)
    folds = splitfolds.(length.(docrng),k);
    #foldrng = mapslices(nd_to_docrng,[countnz(folds[i].==j) for i=1:length(docrng), j=1:k],1);
    return folds
end
