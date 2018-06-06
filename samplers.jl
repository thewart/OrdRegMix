function sample_β{T<:Real}(y,X,σ2_β::T,σ2)
    p = size(X,2);
    Σ = Diagonal(fill(inv(σ2_β),p));
    return sample_β(y,X,Σ,σ2);
end

function sample_β{T<:AbstractVector}(y,X,σ2_β::T,σ2)
    Σ = Diagonal(σ2_β);
    return sample_β(y,X,Σ,σ2);
end

function sample_β{T<:AbstractMatrix}(y,X,Σ_β::T,σ2)
    p = size(X,2);
    Lβ = inv(chol(X'X./σ2 + Σ_β));
    return Lβ*Lβ'*X'y./σ2 + Lβ*randn(p);
 end

function sample_u!(u,σ2_u,y,Z,σ2)
    v = length(Z);
    if v>1
        for i in 1:v
            otherz = setdiff(1:v,i);
            resid = y - Ztu(Z[otherz],u[otherz]);
            sample_u_inner!(u[i],σ2_u[i,:],resid,Z[i],σ2);
        end
    else
        sample_u_inner!(u[1],σ2_u[1,:],y,Z[1],σ2);
    end
end

function sample_u_inner!(u,σ2_u,y,Z,σ2)
    n = size(Z,2);
    d = size(y,2);
    for i in 1:n
        tmp = y[Z[:,i],:];
        X = repeat(eye(d),inner=(size(tmp,1),1));
        u[i,:] .= sample_β(vec(tmp),X,σ2_u.*σ2,σ2);
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

function sample_z(η,y,δ2,cp)
    γ = [-Inf,0,cp,Inf];
    return rand( Truncated( Normal(η,sqrt(δ2)), γ[y],γ[y+1]) );
end

function Ztu{T<:AbstractMatrix}(Z,u::Vector{T})
    return reduce(+,Z.*u)
end
