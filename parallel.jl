using Distributed

function data_step(n::Int64,folder::String,
                   L::Float64,R::Float64,h::Float64,gap::Float64,mu::Float64,p::Float64,con::Float64)
    print("$n)")
    nswaps, Es_opt, dops = simulation(L,R,h,gap,mu,p,con)
    pos = positions(dops)
    occ = [dops[i].occ for i in eachindex(dops)]
    typ = [dops[i].typ for i in eachindex(dops)]
    writedlm(string("/home/usoltcev/Documents/",folder,"/energies/Es$n"),Es_opt)
    writedlm(string("/home/usoltcev/Documents/",folder,"/dopantPos/pos$n"),pos)
    writedlm(string("/home/usoltcev/Documents/",folder,"/dopantOcc/occ$n"),occ)
    writedlm(string("/home/usoltcev/Documents/",folder,"/dopantTyp/typ$n"),typ)
    return nothing
end

function run_data_step(params::Tuple{Int64,String,Float64,Float64,Float64,Float64,Float64,Float64,Float64})
    n,folder,L,R,h,gap,mu,p,con = params
    @time data_step(n,folder,L,R,h,gap,mu,p,con)
    return nothing
end

function data_parallel(set::Int64,folder::String,
                       L::Float64,R::Float64,h::Float64,gap::Float64,mu::Float64,p::Float64,con::Float64)
    ns = 1:set
    params = Iterators.product(ns,[folder],L,R,h,gap,mu,p,con)
    pmap(run_data_step, params)
    
    println("$set simulations completed, please have your files :)")    
    return 1
end