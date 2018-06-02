#read in data
path = "/home/seth/analysis/OrdRegMix/";
v = sum(getindex.(readdir(path),[1:2]) .== "Xr");
Xf = readcsv(string(path,"Xf.csv"));
Xr = [readcsv(string(path,"Xr",i,".csv")) for i=1:v];
Y = readcsv(string(path,"Y.csv"),Int64);
nd = readcsv(string(path,"docrng.csv"),Int64)[:];
docrng = nd_to_docrng(nd);


#read in initialization
init = HYBRIDsample();
#init = init_params(size(Xf,2),size.(Xr,2),size(Y,2));
init.β = readcsv(string(path,"fixef_0.csv"));
init.σ2 = vec(readcsv(string(path,"sigma_0.csv")).^2);
uinit = [readcsv(string(path,"ranef",i,"_0.csv")) for i=1:v];
init.u = [uinit[i][:,d] for i=1:v,d=1:size(Y,2)];
init.α = readcsv(string(path,"alpha_0.csv"));
#init.σ2_u = var.(init.u./sqrt.(init.σ2'));
init.σ2_u = var.(init.u);
init.σ2_u[init.σ2_u.<0.01] = 0.01;

K = size(init.α,2);

using Distributions
resp = readcsv(string(path,"r_0.csv"));
initin = init_params(size(init.α,2),length(docrng),1);
initin.z = vec(rand.(mapslices(Categorical,resp,2)));


#Xf = zeros(size(Xf,1),1);
#Xr = [zeros(size(Xf,1),2)];
Xin = zeros(length(docrng),1);
#hy = hyperparameter(ν0_τ=10.0,τ0_τ=1.0,τ_β=1e-12,ν0_u=1e12,τ0_u=1e-12);
@time lmmtopic(Y,Xf,Xr,Xin,docrng,init=init,initin=initin,K,iter=10,thin=1);
foo = lmmtopic(Y,Xf,Xr,Xin,docrng,init=init,initin=initin,K,iter=500,thin=1);

init.σ2[:] = 1.0;
cp = inv.(vec(readcsv(string(path,"sigma_0.csv"))));
@time foof = lmmtopic_fixedvar(Y,Xf,Xr,Xin,docrng,init=init,
    initin=initin,K,iter=5000,thin=5,cp=cp);
lp = lppd(Y,Xf,Xr,docrng,cp,foof);

initin = init_params(size(init.α,2),1,1);
initin.z = vec(rand.(mapslices(Categorical,resp,2)));
@time foof0 = lmmtopic_fixedvar(Y,Xf,Xr,zeros(1,1),[1:size(Y,1)],init=init,
        initin=initin,K,iter=5000,thin=5,cp=cp);
lp0 = lppd(Y,Xf,Xr,docrng,cp,foof0);

pii = mapslices.(softmax,[foof[i].tlmm.η for i=1:iter],1);
mean(cat(3,pii[501:iter]...),3)
