mutable struct HYBRIDsample
    α::Matrix{Float64}
    β::Matrix{Float64}
    u::Matrix{Vector{Float64}}

    σ2::Vector{Float64}
    σ2_u::Matrix{Float64}

    tlmm::TLMMsample

    HYBRIDsample() = new();
end

function init_params!(samp::HYBRIDsample,p::Int,l::Vector{Int},dim::Int)
    samp.β = randn(p,dim);
    samp.σ2 = rand(dim);
    samp.u = [randn(i) for i=l, d=1:dim];
    samp.σ2_u = rand(v,dim)*2;

    return samp
end

function init_params(p::Int,l::Vector{Int},dim::Int)
    return init_params!(HYBRIDsample(),p,l,dim)
end

function gf(value::Vector{HYBRIDsample},name::Symbol)
    nd = ndims(getfield(value[1],name));
    return cat(nd+1,getfield.(value,name)...)
end

function lmmtopic(y,Xf,Xr,Xin,docrng,K,hy=hyperparameter();hyin=hyperparameter(),iter=1000,thin=1)

    p = size(Xf,2);
    n,dim = size(y);
    v = length(Xr);
    l = size.(Xr,2);

    s = init_params(p,l,dim);
    z = [sample_z(0.0,y[i,d],s.σ2[d]) for i=1:n, d=1:dim];
    fit = TLMMfit([init_params(K,length(docrng),size(Xin,2))],
    Matrix{Vector{NormalMeanPosterior}}(0,0),
    [NormalMeanPosterior(0.0,0.001,s.σ2[d]) for d=1:dim],
    hyin);
    resid = Array{Float64,2}(nobs,dim);

    saveiter = thin:thin:iter;
    nsave = length(saveiter);
    iter = maximum(saveiter);

    samples = Vector{HYBRIDsample}(nsave);

    ##### run model
    for t in 1:iter

        for d in 1:dim
            resid[:,d] = z[:,d] - Xf*s.β[:,d] - Ztu(Xr,s.u[:,d]);
            fit.prior[d].σ2 = s.σ2[d];
        end
        s.α, fit = sample_α(resid',Xin',fit,docrng,K);
        α_expand =  s.α[:,fit.θ[1].z]';

        for d in 1:dim
            resid[:,d] = z[:,d] - (α_expand[:,d] + Ztu(Xr,s.u[:,d]));
            s.β[:,d] = sample_β(resid[:,d],Xf,hy[:τ_β]*s.σ2[d],s.σ2[d]);

            α_Xβ = α_expand[:,d] + Xf*s.β[:,d];
            resid[:,d] = z[:,d] - α_Xβ;
            sample_u!(s.u[:,d],s.σ2_u[:,d].*s.σ2[d],resid[:,d],Xr,s.σ2[d]);
            s.σ2_u[:,d] = sample_σ2.(s.u[:,d]./sqrt(s.σ2[d]),hy[:ν0_u],hy[:τ0_u]);

            α_Xβ_Zu = α_Xβ + Ztu(Xr,s.u[:,d]);
            z[:,d] = sample_z.(α_Xβ_Zu,y[:,d],s.σ2[d]);

            resid[:,d] = z[:,d] - α_Xβ_Zu;
            bigvec = vcat(resid[:,d], s.β[:,d]./sqrt(hy[:τ_β]),
                vcat((s.u[:,d]./sqrt.(s.σ2_u))...));
            s.σ2[d] = sample_σ2(bigvec,hy[:ν0_σ2η],hy[:τ0_σ2η]);

        end

        if t ∈ saveiter
            j = findin(saveiter,t)[1];
            samples[j] = deepcopy(s);
            samples[j].tlmm = deepcopy(fit.θ[1]);
        end
    end

    return samples

end

###### prediction for new samples from old groups
function lppd{U<:PostPredSS}(y,Xf,Xr,fit::Vector{HYBRIDsample{U}})
    iter = length(fit);
    dim = size(y,2);
    lp = zeros(iter);

    for i in 1:iter
        resid .= y .- fit[i].α[:,fit[i].tlmm.z]' .- Xf*fit[i].β[:,d];
        for d in 1:dim
            resid[:,d] .-= Ztu(Xr,fit[i].u[:,d]);
            lp[i] += sum(-0.5 * (resid[:,d]./fit[i].σ2[d]).^2 .- sqrt(fit[i].σ2[d]));
        end
    end

end

function splitfolds(n,k)
    nrep = ceil(Int,n/k);
    return shuffle!(repeat(1:k,inner=nrep))[1:n]
end

function foldup(docrng,k)
    folds = splitfolds.(length.(docrng),k);
    flatfolds = vcat(folds...);
    foldrng = mapslices(nd_to_docrng,[countnz(juh[i].==j) for i=1:length(docrng), j=1:k],1);
    return flatfolds, foldrng
end
