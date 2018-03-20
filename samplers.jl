function sample_β(y,X,σ2_β,σ2)
    p = size(X,2);
    Lβ = inv( chol(X'X./σ2 + I*inv(σ2_β)) );
    return Lβ*Lβ'*X'y./σ2 + Lβ*randn(p);
 end

function sample_u!{T<:AbstractArray}(u::Vector{T},σ2_u,y,X,σ2)
    v = length(X);
    if v>1
        for i in 1:v
            otherz = setdiff(1:v,i);
            u[i] = sample_β(y-Ztu(X[otherz],u[otherz]),X[i],σ2_u[i],σ2);
        end
    else
        u[1] = sample_β(y,X[1],σ2_u[1],σ2);
    end
end

function sample_σ2(y,ν0,τ0);
  a = 0.5*(length(y)+ν0);
  b = 0.5*(τ0*ν0 + dot(y,y));
  return rand(InverseGamma(a,b));
end

function sample_α(y,X,fit,docrng,K)
    fit = topiclmm(y,X,docrng,fit.prior,K,init=fit.θ[1],hy=fit.hyperparameter,iter=1);
    topic = fit.tss[:,:,1];
    α = hcat(rand.(topicpd.(topic))...);
    return α, fit
end

function nd_to_docrng(nd)
    bounds = vcat(1,cumsum(nd));
    n = length(nd);
    return [bounds[i]:bounds[i+1] for i=1:n];
end

function Ztu{T<:AbstractVector}(Z,u::Vector{T})
    return reduce(+,Z.*u)
end

function Ztu{T<:AbstractVector}(Z,u::T)
    return Z*u
end
