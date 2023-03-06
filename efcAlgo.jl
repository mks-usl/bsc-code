# function to check the Efros Shklovskii criterion
function efc(Ej::Float64,Ei::Float64,Vij::Float64)
    efc = false
    if Ej - Ei - Vij > 0
        efc = true
    end
    return efc
end

function swap_occ!(dop1::Dopant,dop2::Dopant)
    occ = dop1.occ
    dop1.occ = dop2.occ
    dop2.occ = occ
    return nothing
end

function swap_counter(con::Int64,result::Int64)
    if result == 1
        con += 1
    end
    return con
end

# return index of random occupied dopant
function randocc(dops::Array{Dopant,1},N::Int64)
    a = 0
    while a == 0
        n = rand(1:N)
        if dops[n].occ != 1
            0
        else 
           return n
        end 
    end
end
# return index of random empty dopant
function randempty(dops::Array{Dopant,1},N::Int64)
    a = 0
    while a == 0
        n = rand(1:N)
        if dops[n].occ != 0
            0
        else 
           return n
        end 
    end
end

# return 0 if no swaps possible, 1 if there is at least one possible
function remain_scr(N::Int64,dops::Array{Dopant,1},Es_scr::Array{Float64,1},Vdops::AbstractArray{Float64,2})
    rem = 0
    for i in 1:N
        for j in 1:N
            if dops[i].occ != 0 || dops[j].occ != 1
            elseif efc(Es_scr[j],Es_scr[i],Vdops[i,j]) == false
            else
                rem = 1
                return rem
            end
        end
    end
    return rem
end

# return number of possible swaps at every given moment
function counter(N::Int64,dops::Array{Dopant,1},Es_scr::Array{Float64,1},Vdops::AbstractArray{Float64,2})
    count = 0
    for i in 1:N
        for j in 1:N
            if dops[i].occ != 0 || dops[j].occ != 1
            elseif efc(Es_scr[j],Es_scr[i],Vdops[i,j]) == false
            else
                count += 1
            end
        end
    end
    
    return count
end

# function trying to optimise a given pair of dopants (i,j) via swapping
function try_to_opt_screened!(i::Int64,j::Int64,
                              dops::Array{Dopant,1},pos::Array{Float64,2},
                              Vdops::Array{Float64,2},Vs_scr::Array{Float64,1},Es_scr::Array{Float64,1},
                              V0_vec::Array{Float64,1},P::Array{Float64,2},
                              zs::Array{Float64,1},phis::Array{Float64,1},
                              gap::Float64,mu::Float64,con::Float64,L::Float64,R::Float64)
    
    # check the Efros Shklovskii criterion first
    if efc(Es_scr[j],Es_scr[i],Vdops[i,j])
            return 0 # no swap needed
        else
        
        if dops[i].occ == dops[j].occ # both filled or empty -> return 0
            return 0
        end
    
        swap_occ!(dops[i],dops[j])

#         Vlin, Slin = linearised(mu,con,V0_vec,P,length(zs),length(phis))
#         Vs_scr = screen_pots(pos,Slin,zs,phis,L,R)
        
        # change energies of the whole configuration
        for k in eachindex(Es_scr)
            if k == i
                Es_scr[k] = dop_E_screened(i,dops,Vdops,Vs_scr,gap)
            elseif k == j
                Es_scr[k] = dop_E_screened(j,dops,Vdops,Vs_scr,gap)
            else
                Es_scr[k] += - Vdops[k,i] + Vdops[k,j]
            end
        end
        
    end
        
    return 1
end

# function to optimise the whole dopant configuration via swapping random pairs
function optimise_scr!(dops::Array{Dopant,1},pos::Array{Float64,2},
                       Vdops::Array{Float64,2},Vs_scr::Array{Float64,1},Es_scr::Array{Float64,1},
                       V0_vec::Array{Float64,1},P::Array{Float64,2},
                       zs::Array{Float64,1},phis::Array{Float64,1},
                       gap::Float64,mu::Float64,con::Float64,L::Float64,R::Float64,p::Float64)
    N = length(dops)
    
    nswaps = 0
    check1 = 0
# as long as there are possible swaps
    while remain_scr(N,dops,Es_scr,Vdops) > 0
        check = 0
        # optimise bulk
        while remain_scr(N,dops,Es_scr,Vdops) > 0
            for i in 1:N^2
                result = try_to_opt_screened!(randocc(dops,N),randempty(dops,N),
                                              dops,pos,Vdops,Vs_scr,Es_scr,
                                              V0_vec,P,zs,phis,gap,mu,con,L,R)
                nswaps = swap_counter(nswaps,result)
            end
        
            check += 1
#             println("$check while-loops completed")
            if check > 10
                break
            end
        end
    
        # update surface
        Es_scr -= Vs_scr
        V, S = compute_surf(mu,p,V0_vec,P,con,length(zs),length(phis))
        Vs_scr = screen_pots(pos,S,zs,phis,L,R)
        Es_scr += Vs_scr
        
#         println("Bulk & surface optimised, check remain...")
	check1 += 1
# safety check to avoid an eternal loop
        if check1 > 10
            break
        end
    end

    return dops,nswaps
end

# function running the full optimisation process after the system initialisation
function run_Efros_screened!(dops::Array{Dopant,1},V0_vec::Array{Float64,1},P::Array{Float64,2},
                             zs::Array{Float64,1},phis::Array{Float64,1},
                             L::Float64,R::Float64,h::Float64,gap::Float64,mu::Float64,con::Float64,p::Float64)
    
    N = length(dops); Nz = length(zs); Nphi = length(phis)
    pos = positions(dops)
    Vdops = pots(dops,L,h)
    
    V, S = compute_surf(mu,p,V0_vec,P,con,Nz,Nphi)
    Vs_scr = screen_pots(pos,S,zs,phis,L,R);
    
    E_scr = Es_screened(dops,Vdops,Vs_scr,gap);
    
    # check for possible swaps
    count = counter(N,dops,E_scr,Vdops) 
    if count == 0
        printstyled("No swaps needed",color=:red); println()
        return 0,E_scr,dops
    end
    
    E_start = total_E_scr(dops,Vdops,Vs_scr,gap)
    
    dops,nswaps = optimise_scr!(dops,pos,Vdops,Vs_scr,E_scr,
                                      V0_vec,P,zs,phis,gap,mu,con,L,R,p)
    
    E_end = total_E_scr(dops,Vdops,Vs_scr,gap)
    println("Energy gain: ", E_start - E_end)
    
    return nswaps,E_scr,dops
end







