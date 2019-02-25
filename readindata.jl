#read in data
v = sum(getindex.(readdir(path),[1:2]) .== "Xr");
Xf = readcsv(string(path,"Xf.csv"));
Xr = [readcsv(string(path,"Xr",i,".csv")) for i=1:v];
Xr = [convert.(Bool,X) for X in Xr];
Y = Bool.(readcsv(string(path,"Y.csv"),Int64));
nd = readcsv(string(path,"docrng.csv"),Int64)[:];
docrng = nd_to_docrng(nd);
Xin = zeros(length(docrng),1);
