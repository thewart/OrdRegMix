function lmmtopic_fixedvar(y,Xf,Xr,Xin,docrng,K;hy=hyperparameter(),hyin=hyperparameter(),
    init=init_params(size(Xf,2),size.(Xr,2),size(y,2)),
    initin=init_params(K,length(docrng),size(Xin,2)),
    cp=repeat(1,outer=(1,size(y,2))),
    iter=1000,thin=1)

    p = size(Xf,2);
    n,dim = size(y);
    v = length(Xr);
    l = size.(Xr,2);

    ν0_α = 0.001;
    μ0_α = 0.0;
    #s = init_params(p,l,dim);
    s = deepcopy(init);

    α_expand =  s.α[:,initin.z]';
    η = Matrix{Float64}(size(y));
    for d=1:dim η[:,d] = α_expand[:,d] + Xf*s.β[:,d] + Ztu(Xr,s.u[:,d]); end
    z = [sample_z(η[i,d],y[i,d],s.σ2[d],cp[d]) for i=1:n, d=1:dim];

    fit = TLMMfit([initin],
        Matrix{Vector{NormalMeanPosterior}}(0,0),
        [NormalMeanPosterior(μ0_α,ν0_α,1.0) for i=1:dim],
        hyin);
    resid = Array{Float64,2}(n,dim);

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
            s.β[:,d] = sample_β(resid[:,d],Xf,hy[:τ_β],s.σ2[d]);

            α_Xβ = α_expand[:,d] + Xf*s.β[:,d];
            resid[:,d] = z[:,d] - α_Xβ;
            sample_u!(s.u[:,d],s.σ2_u[:,d],resid[:,d],Xr,s.σ2[d]);
            s.σ2_u[:,d] = sample_σ2.(s.u[:,d],hy[:ν0_u],hy[:τ0_u]);

            α_Xβ_Zu = α_Xβ + Ztu(Xr,s.u[:,d]);
            z[:,d] = sample_z.(α_Xβ_Zu,y[:,d],s.σ2[d],cp[d]);

            resid[:,d] = z[:,d] - α_Xβ_Zu;
            bigvec = vcat(resid[:,d],(s.α[d,:].-μ0_α)*sqrt(ν0_α),
                    s.β[:,d]./sqrt(hy[:τ_β]),
                    vcat(s.u[:,d]./sqrt.(s.σ2_u[d])...));
            #s.σ2[d] = sample_σ2(bigvec,hy[:ν0_τ],hy[:τ0_τ]);

        end

        if t ∈ saveiter
            j = findin(saveiter,t)[1];
            samples[j] = deepcopy(s);
            samples[j].tlmm = deepcopy(fit.θ[1]);
        end
    end

    return samples

end
