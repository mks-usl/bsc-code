# create type "Dopant"
mutable struct Dopant
    pos::Array{Float64,1}
    typ::Int
    occ::Int
end

# change the representation of an object of type "Dopant"
Base.show(io::IO, dop::Dopant) = print(io, "Dopant at (ρ,ϕ,z) = ", "(", dop.pos[1],",", dop.pos[2],",", dop.pos[2], ")",
                                            " of type ", dop.typ, ", occupied = ", dop.occ )

# create a randomly placed dopant of certain type and occupation
function set_dop(typ::Int64,occ::Int64,L::Float64,R::Float64,h::Float64)
    pos = [sqrt(rand(Uniform(0,1))) * R, (rand(Uniform(0,2π))), rand(Uniform(-L/2, L/2))]
    dop = Dopant(pos,typ,occ)
    return dop
end

# set dopant type
function set_typ(i::Int64)
    if isodd(i)
        typ = 1 # donor
    else
        typ = -1 # acceptor
    end
    return typ
end

# set dopant occupation
function set_occ(i::Int64)
    if isodd(i)
        occ = 0
    else
        occ = 1
    end
    return occ
end

# create a random dopant configuration (charge neutral)
function seed_dops(L::Float64,R::Float64)
    N = 2 * floor(Int, π * R^2 * L)
    
    dops = [set_dop(set_typ(i), set_occ(i), L, R, h) for i in 1:N]
    
    return N,dops
end

# compute dopant charge
function dop_charge(dop::Dopant)
    return ((dop.typ + 1) / 2 - dop.occ) * e
end

# distance between 2 dopants
function dop_dist(dop1::Dopant,dop2::Dopant,L::Float64,h::Float64)
    ri = dop1.pos; rj = dop2.pos
    d = polar_dist(ri,rj,L,h)
    return d
end

# return dopant positions from a dopant array
function positions(dops::Array{Dopant,1})
    N = length(dops)
    pos = zeros(3,N)
    for i in 1:N
        pos[:,i] = dops[i].pos
    end
    return pos
end

# a dopant SPE including surface screening
function dop_E_screened(i::Int64,dops::Array{Dopant,1},Vs::AbstractArray{Float64,2},screenV::Array{Float64,1},gap::Float64)
    sE = 0.5 * gap * dops[i].typ - sum(Vs[i,j] * dop_charge(dops[j]) for j in skip(i,length(dops)))
    return sE - screenV[i]
end

# compute all screened SPEs
function Es_screened(dops::Array{Dopant,1},Vs::AbstractArray{Float64,2},screenV::Array{Float64,1},gap::Float64)
    Es = zeros(length(dops))
    for i in eachindex(dops)
        Es[i] = dop_E_screened(i,dops,Vs,screenV,gap)
    end
    return Es
end

# a technical function to skip realise the summation condition "i != j"
function skip(i,N)
    js = collect(1:N)
    filter!(js -> js != i, js)
    return js
end

# Coulomb potential between two dopants
function V_ij(i::Int64,j::Int64,dops::Array{Dopant,1},L::Float64,h::Float64)
    ab = 1. # eff. Bohr radius
    d = sqrt(dop_dist(dops[i],dops[j],L,h)^2 + ab^2)
    return 1 / d 
end

# compute all potentials between dopants
function pots(dops::Array{Dopant,1},L::Float64,h::Float64)
    N = length(dops)
    Vs = zeros(N,N)
    
    for j in 1:N
        for i in 1:j
            Vs[i,j] = V_ij(i,j,dops,L,h)
        end
    end
    
    return collect(Symmetric(Vs))
end

# compute total energy including screening (without surface kinetic energy)  (optional function)
function total_E_scr(dops::Array{Dopant,1},Vdops::AbstractArray{Float64,2},Vscreen::Array{Float64,1},gap::Float64)
    N = length(dops)
    
    H1 = gap * sum( dops[i].typ * dops[i].occ for i in 1:N )
    H2 = 0.
    for i in 1:N
        for j in skip(i,N)
            H2 += Vdops[i,j] * dop_charge(dops[i]) * dop_charge(dops[j])
        end
    end
    H3 = sum(dop_charge(dops[i]) * Vscreen[i] for i in 1:N) / 2
            
    return (H1 + H2) / 2 + H3
end


