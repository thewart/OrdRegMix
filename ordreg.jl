import DataArrays.cut
using FreqTables.freqtable
using DataFrames, Distributions
n = 500;
p = 4;
α0 = -0.5;
β0 = randn(p);
Xf = randn(p,n);
cp = [-Inf, 0, 2.0, Inf];
Xr = [ModelMatrix(ModelFrame(@formula(y~0+ X),DataFrame(y=fill(0,n),
  X=@pdata(repeat(vcat(1:50),inner=div(n,50)))))).m'];
σ2_u = 1;
u0 = sqrt(σ2_u) * randn(50);
v = length(Xr);

y = cut(α0 + Xf'β0 + Xr[1]'u0 + randn(n),cp);
setlevels!(y,map(string,1:(length(cp)-1)));
y = map(parse,y);

iter = 1000;
prior[:σ2_β] = 5.0;
prior[:σ2_α] = 5.0;
prior[:τ0_δ] = 0.01;
prior[:ν0_δ] = 0.01;
prior[:τ0_u] = 0.01;
prior[:ν0_u] = 0.01;
γ = [-Inf,0,1,Inf];
l = size.(Z,1);

samp = Dict{Symbol,AbstractArray}();
samp[:β] = Array{Float64,3}(p,dim,iter);
samp[:α] = Array{Float64,2}(dim,iter);
samp[:δ2] = Array{Float64,2}(dim,iter);
samp[:u] = Array{Vector{Float64},3}(v,dim,iter);
samp[:σ2_u] = Array{Float64,3}(v,dim,iter);

α = randn();
β = randn(p);
δ2 = 1.0;
u = [randn(l[1])];
σ2_u = rand(v)*2;
z = Vector{Float64}(n);
for i in 1:n
  z[i] = rand( Truncated( Normal(0,2), cp[y[i]],cp[y[i]+1] ));
end

for t in 1:iter
  for d in 1:dim
    z[:,d], α[d], β[:,d], u[:,d], σ2_u[:,d], δ2[d] =
          sample_rate(y[d],z[:,d],Xf,Xr,α[d],β[:,d],u[:,d],σ2_u[:,d],γ,δ2[d],prior)
  end

  samp[:β][:,:,t] = β;
  samp[:α][:,t] = α;
  samp[:δ2][:,t] = δ2;
  samp[:u][:,:,t] = deepcopy(u);
  samp[:σ2_u][:,:,t] = σ2_u;
end


function calcZtu(Z,u)
  if length(Z)>0
    Ztu = fill(0.0,size(Z[1],2));
    for i in 1:length(Z) Ztu += Z[i]'u[i]; end
  else
    Ztu=0;
  end
  return Ztu
end

function sample_rate(y,z,Xf,Xr,α,β,u,σ2_u,γ,δ2,prior)

  α = sum(z-Xf'β-calcZtu(Xr,u))/(n+inv(prior[:σ2_α])) +
      sqrt(δ2/(n+inv(prior[:σ2_α]))) * randn();

  Lβ = inv( chol(Xf*Xf' + I*inv(prior[:σ2_β])) );
  β = Lβ*lβ'*Xf*(z.-α-calcZtu(Xr,u)) + sqrt(δ2)*Lβ*randn(p);

  for i in 1:v
    otherz = setdiff(1:v,i);
    Lu = inv( chol(Xr[i]*Xr[i]' + I*inv(σ2_u[i])) );
    u[i] = Lu*Lu'*Xr[i]*(z.-α-Xf'β-calcZtu(Z[otherz],u[otherz])) +
     sqrt(δ2)*Lu*randn(l[i]);

    a = 0.5(l[i]+prior[:ν0_u]);
    b = 0.5(prior[:σ0_u]*prior[:ν0_u] + dot(u[i],u[i])/δ2);
    σ2_u[i] = rand(InverseGamma(a,b));
  end

  η = α + Xf'β + calcZtu(Xr,u);
  for i in 1:n
    z[i] = rand( Truncated( Normal(η[i],sqrt(δ2)), γ[y[i]],γ[y[i]+1]) );
  end

  ss = z - α - Xf'β - calcZtu(Xr,u);
  a = 0.5*(n+prior[:ν0_δ]);
  b = 0.5*(prior[:δ0_δ]*prior[:ν0_δ] + dot(ss,ss));
  δ2 = rand(InverseGamma(a,b));

  return z, α, β, u, σ2_u, δ2
end
