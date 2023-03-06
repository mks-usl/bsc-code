using Distributed

hosts = [("l85",4)]
addprocs(hosts, topology=:master_worker, exename=`nice -n 19  $(Sys.BINDIR)/julia`)
@everywhere begin
	include("include.jl")
	include("systemParams.jl")
end
folder = "SimParallelR7L100/gap10/Mu0/Simulation7test";
set = 60
@time data_parallel(set,folder,L,R,h,gap,mu,p,con)




