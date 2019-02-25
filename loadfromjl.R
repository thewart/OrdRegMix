library(JuliaCall)
source('~/code/MultVarBinom/helperfuncs.R')

if (!exists("julia")) {
  julia <- julia_setup(JULIA_HOME = "/home/seth/julia/usr/bin/")
  julia_command("include(\"/home/seth/code/OrdRegMix/initialize.jl\")")
}
julia_assign("path",paste0("/home/seth/analysis/OrdRegMix/",path))
julia_command("include(\"/home/seth/code/OrdRegMix/readindata.jl\");")

julia_command("@load string(path,\"topicfits.jld\") foof foofu foof0 lp lpu lp0;")
