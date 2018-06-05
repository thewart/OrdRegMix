function sample_β(y::Vector{Float64},X,σ2_β,σ2)
    p = size(X,2);
    Lβ = inv( chol(X'X./σ2 + I*inv(σ2_β)) );
    return Lβ*Lβ'*X'y./σ2 + Lβ*randn(p);
 end

function sample_u_inner(y,Λ,σ2)
    n,d = size(y);
    Lu = inv(chol(Diagonal(inv.(σ2).*n) + Λ));
    return Lu*Lu'*vec(sum(y,1))./σ2 + Lu*randn(d);
end

function sample_u!(u,y,Z,Σ,σ2)
    Λ = Hermitian(inv(Σ));
    for i in 1:size(Z,2)
        u[i,:] = sample_u_inner(y[Z[:,i],:],Λ,σ2);
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

function sample_Σ(u,σ2diag,ν0,a0)
    a = rand.(InverseGamma.(0.5*(ν0+length(σ2diag)),ν0.*σ2diag .+ a0^-2))
    ν = ν0 + size(u,1) + size(u,2);
    P = u'u + Diagonal(inv.(a));
    return rand(InverseWishart(ν,P))
end

function Ztu{T<:AbstractArray}(Z,u::Vector{T})
    return reduce(+,Z.*u)
end

function Ztu{T<:AbstractMatrix}(Z,u::T)
    return Z*u
end
