import DataArrays.cut
using FreqTables.freqtable
using DataFrames, Distributions
n = 500;
p = 4;
α0 = -0.5;
β0 = randn(p);
X = randn(p,n);
cp = [-Inf, 0, 2.0, Inf];

y = cut(α0 + X'β0 + randn(n),cp);
setlevels!(y,map(string,1:(length(cp)-1)));
y = map(parse,y);

iter = 1000;
τ_β = 5.0;
ν0 = 0.01;
δ0 = 0.01;
γ = [-Inf,0,1,Inf];
Lβ = inv( chol( X*X' + I*inv(τ_β)));
ΣβX = Lβ*Lβ'*X;
samp = Dict{Symbol,AbstractArray}();
samp[:β] = Array{Float64,2}(p,iter);
samp[:α] = Vector{Float64}(iter);
samp[:δ2] = Vector{Float64}(iter);

α = randn();
β = randn(p);
δ2 = 1.0;
z = Vector{Float64}(n);
for i in 1:n z[i] = rand( Truncated( Normal(0,2), cp[y[i]],cp[y[i]+1] ));
end

for t in 1:iter
  α = sum(z.-X'β)/(n+ν0) + sqrt(δ2)/(n+ν0) * randn();

  β = ΣβX*(z.-α) + sqrt(δ2)*Lβ*randn(p);

  η = α + X'β;
  for i in 1:n z[i] = rand( Truncated( Normal(η[i],sqrt(δ2)), γ[y[i]],γ[y[i]+1]) );
  end

  ss = z - α - X'β;
  a = 0.5*(n+ν0);
  b = 0.5*(δ0*ν0 + dot(ss,ss));
  δ2 = rand(InverseGamma(a,b));

  samp[:β][:,t] = β;
  samp[:α][t] = α;
  samp[:δ2][t] = δ2;
end
