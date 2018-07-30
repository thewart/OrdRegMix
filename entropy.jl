function calc_cond_entropy(fit::HYBRIDsample,n::Int=200)
    dim = size(fit.α,1);
    K = size(fit.α,2);
    L = chol(fit.σ2_u[1])';
    cH_pi = 0.0;
    cH = 0.0;
    epy_pi = fill(0.0,2^dim);
    epy = fill(0.0,2^dim);

    pk = [sample_pk(fit.tlmm,K) for i in 1:n];
    αc = [fit.α .- fit.α*i for i in pk];

    for j in 1:n
        fill!(epy_pi,0.0);
        for i in 1:n
            pyd = calc_pyd(αc[j],L*randn(dim)+fit.α*pk[j],K,dim);
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
    #sample from p(π)
    for i=1:K ηk[i] = tlmm.μ[i] + sqrt(tlmm.σ2_η[i])*randn(); end
    return softmax(ηk)
end

function calc_pyd(α,μ,K,dim)
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
