function calc_entropy_parallel(fit::HYBRIDsample,cp,n::Int)
    dim = length(fit.σ2);
    K = size(fit.α,2);
    L = chol(fit.σ2_u[1])';
    epy = fill(0.0,3^dim);
    cH = 0.0;
    py = SharedVector{Float64}(3^dim);

    for i in 1:n
        pk = sample_pk(fit.tlmm,K);
        pyd = calc_pyd(fit.α,L*randn(dim),K,dim);
        calc_py_parallel!(py,pyd,pk);
        epy .+= py;
        cH += entropy(py);
    end
    cH *= 1/n;
    epy .*= 1/n;
    H = entropy(epy);

    return H, cH
end

function calc_cond_entropy_parallel(fit::HYBRIDsample,n::Int,nin::Int,uin::Bool=false)
    dim = length(fit.σ2);
    K = size(fit.α,2);
    L = chol(fit.σ2_u[1])';
    cH = 0.0;
    py = SharedVector{Float64}(3^dim);
    epy_in = fill(0.0,3^dim);
    epy_out = fill(0.0,3^dim);

    for i in 1:n
        fill!(epy_in,0.0);
        if uin
            pyd = calc_pyd(fit.α,L*randn(dim),fit.σ2,K,dim);
        else
            pk = sample_pk(fit.tlmm,K);
        end
        for j in 1:nin
            if uin
                pk = sample_pk(fit.tlmm,K);
            else
                pyd = calc_pyd(fit.α,L*randn(dim),K,dim);
            end
            calc_py_parallel!(py,pyd,pk);
            epy_in .+= py;
        end
        epy_out .+= epy_in;
        cH += entropy(epy_in./nin);
    end
    cH *= 1/n;
    epy_out .*= 1/(n*nin);
    H = entropy(epy_out);

    return H, cH
end

function calc_cond_entropy_parallel_2(fit::HYBRIDsample,n::Int)
    dim = length(fit.σ2);
    K = size(fit.α,2);
    L = chol(fit.σ2_u[1])';
    cH_pi = 0.0;
    cH = 0.0;
    py = SharedVector{Float64}(3^dim);
    # epy_u = fill(0.0,3^dim,n);
    epy_pi = fill(0.0,3^dim);
    epy = fill(0.0,3^dim);

    pk = [sample_pk(fit.tlmm,K) for i in 1:n];
    αc = [fit.α .- fit.α*i for i in pk];
    # u = [L*randn(dim) for i in 1:n];

    for j in 1:n
        fill!(epy_pi,0.0);
        for i in 1:n
            pyd = calc_pyd(αc[j],L*randn(dim)+fit.α*pk[j],fit.σ2,K,dim);
            calc_py_parallel!(py,pyd,pk[j]);
            cH += entropy(py);
            epy_pi .+= py;
            # epy_u[:,i] .+= py;
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

function calc_cond_entropy_2(fit::HYBRIDsample,n::Int=200)
    dim = length(fit.σ2);
    K = size(fit.α,2);
    L = chol(fit.σ2_u[1])';
    cH_pi = 0.0;
    cH = 0.0;
    epy_pi = fill(0.0,3^dim);
    epy = fill(0.0,3^dim);

    pk = [sample_pk(fit.tlmm,K) for i in 1:n];
    αc = [fit.α .- fit.α*i for i in pk];

    for j in 1:n
        fill!(epy_pi,0.0);
        for i in 1:n
            pyd = calc_pyd(αc[j],L*randn(dim)+fit.α*pk[j],fit.σ2,K,dim);
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

function calc_py_parallel!(py,pyd,pk)
    @sync for j in procs(py)
        @async remotecall_wait(increment_py!,j,py,pyd,pk);
    end
end

function increment_py!(py::SharedVector{Float64},pyd,pk)
    K,dim = size(pyd,2,3);
    tmp = Vector{Float64}(K);

    #compute probabilities of each possible y (oy!)
    for i in localindexes(py)
        for k in 1:K tmp[k] = pk[k]; end
        ind = i-1;
        for d in 1:dim
            ind, rem = divrem(ind,3);
            for k in 1:K
                tmp[k] *= pyd[rem+1,k,d];
            end
        end
        py[i] = sum(tmp);
    end
end

function calc_py(pyd,pk)
    K,dim = size(pyd,2,3);
    tmp = Vector{Float64}(K);
    py = Vector{Float64}(3^dim);
    #compute probabilities of each possible y (oy!)
    for i in 1:length(py)
        for k in 1:K tmp[k] = pk[k]; end
        ind = i-1;
        for d in 1:dim
            ind, rem = divrem(ind,3);
            for k in 1:K
                tmp[k] *= pyd[rem+1,k,d];
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

function calc_pyd(α,μ,σ2,K,dim)
    pyd = Array{Float64,3}(3,K,dim);
    for d in 1:dim
        N = Normal(0,sqrt(σ2[d]));
        for k in 1:K
            η = μ[d] + α[d,k];
            #compute per-behavior probabilities for current η
            pyd[1,k,d] = cdf(N,-η);
            pyd[3,k,d] = cdf(N,η-1.0);
            pyd[2,k,d] = 1 - pyd[1,k,d] - pyd[3,k,d];
        end
    end
    return pyd
end
