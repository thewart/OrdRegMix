function calc_entropy(fit::HYBRIDsample,cp,n::Int)
    py = calc_py(fit,cp,n::Int);
    return -sum(py.*log.(py))
end

function calc_py(fit::HYBRIDsample,cp,n::Int)
    dim = length(fit.σ2);
    K = size(fit.α,2);
    py = zeros(3^dim);

    for i in 1:n
        pk = sample_pk(fit.tlmm,K);
        pyd = calc_pyd(fit,K,dim);

        increment_py!(py,pyd,pk);
    end
    py .*= 1/n;

    return py
end

function sample_pk(tlmm,K)
    ηk = Vector{Float64}(K);
    #sample from p(π)
    for i=1:K ηk[i] = fit.tlmm.μ[i] + sqrt(fit.tlmm.σ2_η[i])*randn(); end
    return softmax(ηk)
end

function calc_pyd(fit,K,dim)
    pyd = Array{Float64,3}(3,K,dim);
    for d in 1:dim
        for k in 1:K
            η = u + fit.α[d,k];
            #compute per-behavior probabilities for current η
            pyd[1,k,d] = normcdf(-η);
            pyd[3,k,d] = normcdf(η-cp[d]);
            pyd[2,k,d] = 1 - pyd[1,k,d] - pyd[3,k,d];
        end
    end
    return pyd
end

function py(pyd,pk)
    K,dim = size(pyd,2,3);
    ny = 3^dim;
    p = Vector{Float64}(ny);
    tmp = Vector{Float64}(K);

    #compute probabilities of each possible y (oy!)
    for i in 1:ny
        for k in 1:K tmp[k] = pk[k]; end
        ind = i-1;
        for d in 1:dim
            ind, rem = divrem(ind,3);
            for k in 1:K
                tmp[k] *= pyd[rem+1,k,d];
            end
        end
        p[i] = sum(tmp);
    end
    return p
end

function increment_py!(py,pyd,pk)
    K,dim = size(pyd,2,3);
    ny = 3^dim;
    tmp = Vector{Float64}(K);

    #compute probabilities of each possible y (oy!)
    for i in 1:ny
        for k in 1:K tmp[k] = pk[k]; end
        ind = i-1;
        for d in 1:dim
            ind, rem = divrem(ind,3);
            for k in 1:K
                tmp[k] *= pyd[rem+1,k,d];
            end
        end
        py[i] += sum(tmp);
    end
end
