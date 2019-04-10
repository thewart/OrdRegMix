@everywhere using DataFrames, Distributions, StatsBase, StatsFuns;
@everywhere include("/home/seth/code/LogisticTopicRegression/LogTopReg.jl");
@everywhere using LogTopReg;
import LogTopReg.init_params, .LogTopReg.init_params!;
@everywhere include("/home/seth/code/OrdRegMix/hybridout.jl");
@everywhere include("/home/seth/code/OrdRegMix/samplers.jl");
include("/home/seth/code/OrdRegMix/entropy.jl");
using JLD;
