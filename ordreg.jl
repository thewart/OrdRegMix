import DataArrays.cut
using FreqTables.freqtable
using DataFrames, Distributions
n = 1000;
p = 4;
v
α0 = -0.25;
β0 = randn(p);
X = randn(p,n);
cp = [-Inf, 0, 2.0, Inf];
Z = [ModelMatrix(ModelFrame(@formula(y~0+ X),DataFrame(y=fill(0,n),
  X=@pdata(repeat(vcat(1:50),inner=div(n,50)))))).m'];
σ2_u = 1.5;
u0 = sqrt(σ2_u) * randn(50);
v = length(Z);

y = cut(α0 + X'β0 + Z[1]'u0 + randn(n),cp);
setlevels!(y,map(string,1:(length(cp)-1)));
y = map(parse,y);

iter = 1000;
τ_β = 5.0;
ν0 = 0.01;
δ0 = 0.01;
σ0_u = 0.01;
ν0_u = 0.01;
γ = [-Inf,0,1,Inf];
l = size.(Z,1);

Lβ = inv( chol(X*X' + I*inv(τ_β)) );
ΣβX = Lβ*Lβ'*X;

samp = Dict{Symbol,AbstractArray}();
samp[:β] = Array{Float64,2}(p,iter);
samp[:α] = Vector{Float64}(iter);
samp[:δ2] = Vector{Float64}(iter);
samp[:u] = Array{Vector{Float64}}(v,iter);
samp[:σ2_u] = Array{Float64,2}(v,iter);

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
  α = sum(z-X'β-calcZtu(Z,u))/(n+ν0) + sqrt(δ2)/(n+ν0) * randn();

  β = ΣβX*(z.-α-calcZtu(Z,u)) + sqrt(δ2)*Lβ*randn(p);

  for i in 1:v
    otherz = setdiff(1:v,i);
    Lu = inv( chol(Z[i]*Z[i]' + I*inv(σ2_u[i])) );
    u[i] = Lu*Lu'*Z[i]*(z.-α-X'β-calcZtu(Z[otherz],u[otherz])) +
     sqrt(δ2)*Lu*randn(l[i]);

    a = 0.5(l[i]+ν0_u);
    b = 0.5(σ0_u*ν0_u + dot(u[i],u[i])/δ2);
    σ2_u[i] = rand(InverseGamma(a,b));
  end

  η = α + X'β + calcZtu(Z,u);
  for i in 1:n
    z[i] = rand( Truncated( Normal(η[i],sqrt(δ2)), γ[y[i]],γ[y[i]+1]) );
  end

  ss = z - α - X'β - calcZtu(Z,u);
  a = 0.5*(n+ν0);
  b = 0.5*(δ0*ν0 + dot(ss,ss));
  δ2 = rand(InverseGamma(a,b));

  samp[:β][:,t] = β;
  samp[:α][t] = α;
  samp[:δ2][t] = δ2;
  samp[:u][:,t] = u;
  samp[:σ2_u][:,t] = σ2_u;
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
