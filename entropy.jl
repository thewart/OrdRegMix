function calc_cond_entropy(fit::HYBRIDsample,n::Int=200,αmean::Bool=false)
    dim = size(fit.α,1);
    K = size(fit.α,2);
    L = chol(fit.σ2_u[1])';
    cH_pi = 0.0;
    cH = 0.0;
    epy_pi = fill(0.0,2^dim);
    epy = fill(0.0,2^dim);

    pk = [sample_pk(fit.tlmm,K) for i in 1:n];
    αc = [fit.α .- fit.α*i for i in pk];
    pμ = softmax(fit.tlmm.μ);

    for j in 1:n
        fill!(epy_pi,0.0);
        for i in 1:n
            pyd = calc_pyd(αc[j],L*randn(dim)+fit.α*pk[j]*!αmean+fit.α*pμ*αmean,K,dim);
            py = calc_py(pyd,pk[j]);
            cH += entropy(py);
            epy_pi .+= py;
        end
        epy .+= epy_pi;
        cH_pi += entropy(epy_pi./n);
    end

    cH_pi *= 1/n;
    cH *=1/n^2;

    # epy_u .*= 1/n;
    # cH_u = mean(mapslices(entropy,epy_u,1));

    epy .*= 1/(n^2);
    H = entropy(epy);

    return H, cH, cH_pi
end

function calc_cond_entropy_2(fit::HYBRIDsample)
    dim = size(fit.α,1);
    K = size(fit.α,2);
    n = size(fit.u[1],1);

    #pk = [sample_pk(fit.tlmm,K) for i in 1:n];
    pi = mapslices(softmax,fit.tlmm.η,1);
    αbar = fit.α * mean(pi,2);
    α = fit.α * (pi .- mean(pi,2));
    u = Eα + fit.u[1]';

    for j in 1:n
        pyd = calc_pyd(αc[j],u[i]+fit.α*pk[j],K,dim);
        py = calc_py(pyd,pk[j]);
            cH += entropy(py);
            epy_col .+= py;
            epy_row[:,i] .+= py;
        end
        epy .+= epy_col;
        cH_col += entropy(epy_col./n);
    end

    cH_col *= 1/n;
    cH *=1/n^2;

    epy_row .*= 1/n;
    cH_row = mean(mapslices(entropy,epy_row,1));

    epy .*= 1/(n^2);
    H = entropy(epy);

    return H
end


function calc_py(pyd,pk)
    K,dim = size(pyd,1,3);
    tmp = Vector{Float64}(K);
    py = Vector{Float64}(2^dim);
    #compute probabilities of each possible y (oy!)
    for i in 1:length(py)
        for k in 1:K tmp[k] = pk[k]; end
        ind = i-1;
        for d in 1:dim
            ind, rem = divrem(ind,2);
            for k in 1:K
                tmp[k] *= pyd[k,rem+1,d];
            end
        end
        py[i] = sum(tmp);
    end
    return py;
end

function sample_pk(tlmm,K)
    ηk = Vector{Float64}(K);
    #estimate effective variance
    η = tlmm.η .- tlmm.η[1:1,:];
    σ2 = mapslices(var,η,2);
    μ = mapslices(mean,η,2);
    #sample from p(π)
    #for i=1:K ηk[i] = tlmm.μ[i] + sqrt(tlmm.σ2_η[i])*randn(); end
    for i=1:K ηk[i] = μ[i] + sqrt(σ2[i])*randn(); end
    return softmax(ηk)
end

function calc_pyd(α,μ,K=size(α,2),dim=size(α,1))
    pyd = Array{Float64,3}(K,2,dim);
    for d in 1:dim
        for k in 1:K
            η = μ[d] + α[d,k];
            #compute per-behavior probabilities for current η
            pyd[k,2,d] = normcdf(η);
            pyd[k,1,d] = 1-pyd[k,2,d];
        end
    end
    return pyd
end

function calc_py(fit::Vector{HYBRIDsample})
    n,d = size(fit[1].u[1]);
    py = Array{Float64}(2^d,n,length(fit));
    for j in eachindex(fit)
        s = fit[j];
        for l in 1:n
            pi = softmax(s.tlmm.η[:,l]);
            u = s.u[1][l,:];
            pyd = calc_pyd(s.α,u);
            py[:,l,j] = calc_py(pyd,pi);
        end
    end
    return py
end
