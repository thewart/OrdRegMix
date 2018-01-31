function sample_β(y,X,σ2_β,δ2)
    Lβ = inv( chol(X'X + I*inv(σ2_β)) );
    return Lβ*Lβ'*X'y + sqrt(δ2)*Lβ*randn(p);
 end

function sample_u!(u,σ2_u,y,X,δ2)
    v = length(X);
    for i in 1:v
        otherz = setdiff(1:v,i);
        u[i] = sample_β(y-Ztu(X[otherz],u[otherz]),X[i],σ2_u[i],δ2);
    end
end

function sample_z(η,y,γ,δ2)
  return rand( Truncated( Normal(η,sqrt(δ2)), γ[y[i]],γ[y[i]+1]) );
end

function sample_σ2(y,ν0,τ0,σ2=1.0);
  a = 0.5*(length(y)+ν0);
  b = 0.5*(τ0*ν0 + dot(y,y)/σ2);
  return rand(InverseGamma(a,b));
end

function sample_α(y,X,fit,topic,K)
    fit, topic = topiclmm(y,X,fit.prior,topic,K,init=fit.θ[1],hy=fit.hyperparameter,iter=1);
    topic = topic[:,:,1];
    α = getfield.(topic,:μ);
    α = α[:,fit.θ[1].z]';
    return α, topic, fit
end
