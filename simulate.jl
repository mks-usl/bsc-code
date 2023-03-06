# initialise the whole system
function initialise(L::Float64,R::Float64,h::Float64)
    # discretise in z and phi direction
    zs = collect(-L/2 + h/2 : h : L/2)
    k = round(2π / h * R); m = 2π / k
    phis = collect(-π + m/2 : m : π)
#     println("Wire initialised with length $L, radius $R and discretisation step $h")
    
    # create randomly distributed dopants
    Ndops = 2*floor(Int, π * R^2 * L)
    dops = [set_dop(set_typ(i), set_occ(i), L, R, h) for i in 1:Ndops]
    qs = dop_charge.(dops)
    pos = positions(dops)
#     println("$Ndops dopants seeded at random positions")
            
    # compute the external potential
    V0 = zeros(length(zs),length(phis))
    for n in eachindex(qs)
        for i in eachindex(zs)
            for j in eachindex(phis)
                V0[i,j] += V_ext(zs[i],phis[j],pos[:,n],qs[n],h,L,R) # pos in polar: (rho,phi,z)
            end
        end
    end
    V0_vec = mat_to_vec(V0)
#     println("External potential for $Ndops charges computed")
    
    # compute the patch integrals on surface
    P = init_patch(zs,phis,h,m,L,R)
#     println("Patch integrals computed")
    
#     printstyled("Initialisation completed, time elapsed:",color=:green)
    return dops,V0_vec,zs,phis,P
end

# run the full simulation starting with the initial parameters
function simulation(L::Float64,R::Float64,h::Float64,gap::Float64,mu::Float64,p::Float64,con::Float64)
    # initialise the wire
    dops,V0_vec,zs,phis,P = initialise(L,R,h)
    Nz = length(zs); Nphi = length(phis)
    V0 = vec_to_mat(V0_vec, Nz, Nphi);
    N = length(dops);
    
    nswaps,Es_opt,dops = run_Efros_screened!(dops,V0_vec,P,zs,phis,L,R,h,gap,mu,con,p);
    
    return nswaps, Es_opt, dops
end

# function to save data into a folder in my home directory
function collect_data(set::Int64,folder::String,
                      L::Float64,R::Float64,h::Float64,gap::Float64,mu::Float64,p::Float64,con::Float64)
    
    for n in 1:set
        print("$n)")
        nswaps, Es_opt, dops = simulation(L,R,h,gap,mu,p,con)
        pos = positions(dops)
        occ = [dops[i].occ for i in eachindex(dops)]
        typ = [dops[i].typ for i in eachindex(dops)]
        writedlm(string("/home/usoltcev/Documents/",folder,"/energies/Es$n"),Es_opt)
        writedlm(string("/home/usoltcev/Documents/",folder,"/dopantPos/pos$n"),pos)
        writedlm(string("/home/usoltcev/Documents/",folder,"/dopantOcc/occ$n"),occ)
        writedlm(string("/home/usoltcev/Documents/",folder,"/dopantTyp/typ$n"),typ)
    end
    
    println("$set simulations completed, please have your files :)")    
    return 1
end












