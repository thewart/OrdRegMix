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
