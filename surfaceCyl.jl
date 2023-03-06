# smooth function for the periodic boundary conditions
f(x,d) = 0.5 * (1 + x / sqrt(x^2 + d^2))
function pbc(z1::Float64,z2::Float64,L::Float64,h::Float64)
    dz = abs(z1 - z2)
    return f(L/2 - dz,2h) * dz + f(-L/2  + dz,2h) * (L - dz)
end

# transform polar of form (rho,phi,z) into cartesian of form (x,y,z)
function p_in_c(r_p::Array{Float64,1})
    ρ,ϕ,z = r_p
    x = ρ * cos(ϕ)
    y = ρ * sin(ϕ)
    return [x,y,z]
end

# distance between two surface points of form (rho,phi,z)
function polar_dist(r::Array{Float64,1},r0::Array{Float64,1},L::Float64,h::Float64)
    ρ,ϕ,z = r ; ρ0,ϕ0,z0 = r0
    x,y = [ρ*cos(ϕ), ρ*sin(ϕ)] ; x0,y0 = [ρ0*cos(ϕ0), ρ0*sin(ϕ0)]
    dz = pbc(z,z0,L,h)
    return sqrt((x-x0)^2 + (y-y0)^2 + dz^2)
end

# transform an array of polar points of form (rho,phi,z) into cartesian positions
function cart_pos(pos::Array{Float64,2})
    N = size(pos,2)
    pos_c = zeros(3,N)
    for i in 1:N
        pos_c[:,i] = p_in_c(pos[:,i])
    end
    return pos_c
end

# fluctuation of the surface charge consentration as function of electrostatic potential and chemical potential
function δn(V::Float64,μ::Float64)
    dn = (V + μ)^2 * sign(V + μ) - μ^2 * sign(μ)
    return dn
end

# transform charge concentration into charge density via multiplying with the electron charge -e
function dens_ij(i::Int64,j::Int64,V::Array{Float64,2},μ::Float64,con::Float64)
    s = - e * con * δn(V[i,j],μ)
    return s
end

# compute the array of surface charge density values at each surface point starting from the electrostatic potential values
function dens(V::Array{Float64,2},μ::Float64,con::Float64)
    
    sigmas = zeros(size(V,1),size(V,2))
    
    for j in 1:size(V,2)
        for i in 1:size(V,1)
            sigmas[i,j] = dens_ij(i,j,V,μ,con)
        end
    end
    
    return sigmas
end

# Coulomb potential between two cartisian points
function V_coul(r::Array{Float64,1},r0::Array{Float64,1},q::Float64,h::Float64,L::Float64)
    ab = 1. # effective Bohr radius
    V_r = q / sqrt( (r[1]-r0[1])^2 + (r[2]-r0[2])^2 + pbc(r[3],r0[3],L,h)^2 + ab )
    return V_r
end

# Coulomb potential (external potential) at a surface point (R,phi,z) for a charge at q_pos in the bulk
function V_ext(z,phi,q_pos::Array{Float64,1},q::Float64,h::Float64,L::Float64,R::Float64)
    r  = p_in_c([R,phi,z])
    r0 = p_in_c(q_pos)
    V = V_coul(r,r0,q,h,L)
    return V
end

# functions to transform a 2d-array into a 1d-list of tuples (simplifies calculation of Coulomb integrals on surface)
function iter_x(i::Int64,N::Int64)
    return div((i-1),N)+1
end

function iter_y(i::Int64,N::Int64)
    if mod(i,N) != 0
        return mod(i,N)
    else
        return N
    end
end

function mat_to_vec(M::Array{Float64,2})
    Nx = size(M,1)
    Ny = size(M,2)
    vec = zeros(Nx*Ny)
    
    for i in eachindex(vec)
        vec[i] = M[iter_x(i,Ny), iter_y(i,Ny)]
    end
    
    return vec
end

function vec_to_mat(vec::Array{Float64,1},Nx::Int64,Ny::Int64) # works only for known Nx, Ny
    M = zeros(Nx,Ny)
    
    for i in 1:Nx
        M[i,:] = vec[Ny*(i-1) + 1 : i*Ny]
    end
    return M
end

# function to check convergence (compare step n with step n+1)
function conv_var(V1_vec::Array{Float64,1},V2_vec::Array{Float64,1})
    s = 0.
    for a in eachindex(V1_vec)
        s += abs(V1_vec[a]^2 - V2_vec[a]^2)
    end
    return sqrt(s) / length(V1_vec)
end

function check_conv(V1_vec::Array{Float64,1},V2_vec::Array{Float64,1},p::Float64)
    var = conv_var(V1_vec,V2_vec)
    converged = false
    if var <= p
        converged = true
    end
    return converged
end

# a Coulomb integral between two points on surface a and b for a list of points
function P_ab(a::Int64,b::Int64,points::Array{Tuple{Float64,Float64},1},h::Float64,m::Float64,R::Float64,l::Float64)
    z,phi = points[b]
    Z,Phi = points[a]
    starts = [Z - h/2, Phi - m/2] 
    ens = [Z + h/2, Phi + m/2]

    patch(r) = @. R / sqrt( pbc(r[1],z,l,h)^2 + 2 * R^2 * (1 - cos(r[2]-phi)) + h^2 * 1e-8)
    P = hcubature(patch,starts,ens,abstol=1e-8)[1]
    
    return P
end

# function using symmetries of system to compute integrals between each and every pair of points from the list
function P_blocks(points::Array{Tuple{Float64,Float64},1},Nphi::Int64,h::Float64,m::Float64,R::Float64,L::Float64)
    N = length(points)
    P = zeros(N,N)
    
    # first Nz blocks
    for b in 1:N
        for a in 1:Nphi
            P[a,b] = P_ab(a,b,points,h,m,R,L)
        end
    end
    
    # fill out the matrix
    for a in Nphi+1:N
        for b in 1:N
            P[a,b] = P[mod(a - (Nphi+1),N)+1, mod(b - (Nphi+1),N)+1]
        end
    end
    
    
    return P
end

# initialise a patch integral matrix from given arrays of z and phi values
function init_patch(zs::Array{Float64,1},phis::Array{Float64,1},h::Float64,m::Float64,L::Float64,R::Float64)
    Nz = length(zs)
    Nphi = length(phis)
    points = [(zs[iter_x(i,Nphi)], phis[iter_y(i,Nphi)]) for i in 1:Nz*Nphi];

    Pb = P_blocks(points,Nphi,h,m,R,L)
    
    return Pb
end

# do the full integral over patches via summing up dn * patch integral between every pair (a,b)
function int_patch(b::Int64,
                   V_vec::Array{Float64,1},
                   P::Array{Float64,2},
                   mu::Float64,con::Float64)
    
    s = 0.
    for a in eachindex(V_vec)
        s += δn(V_vec[a],mu) * P[a,b]
    end
    
    return con * s
end

# one step of the relaxation algorithm
function step_patch(V_vec::Array{Float64,1},V0_vec::Array{Float64,1},
                    P::Array{Float64,2},
                    γ::Float64,mu::Float64,con::Float64
                    )

    V_new = zeros(length(V_vec))

    # compute new potential array
    for b in eachindex(V_new)
        V_new[b] = V0_vec[b] - int_patch(b,V_vec,P,mu,con)
    end
    
    # try to use an optimal gamma (optional)
    g = γ ^ 2
    if (1-g) * maximum(abs.(V_new)) - g * maximum(abs.(V_vec)) < 0
        γ = g
    end
    
    # relaxation step
    V_vec = @. (1-γ) * V_new + γ * V_vec
    
    return V_vec
end

# adjust a value of the relaxation parameter gamma with an empiric fit function
func1(x) = 0.9496417 + (0.2377924 - 0.9496417)/(1 + (x/32.27631)^2.027644)
function g_new(mu::Float64)
    if abs.(mu) >= 10
        g = func1(abs.(mu))
    else
        g = 0.25
    end
    g
end

# do the full relaxation until convergence level reached (precision p)
function conv_patch(p::Float64,V0_vec::Array{Float64,1},P::Array{Float64,2},mu::Float64,con::Float64)
    N_max = 1000
    γ = g_new(mu)

    V_vec = copy(V0_vec)
    
    n = 0
    converged = false
    for i in 1:N_max
        V_before = copy(V_vec)
        V_vec = step_patch(V_vec,V0_vec,P,γ,mu,con)
        
        if mod(i,1) == 0
            converged = check_conv(V_vec,V_before,p)
            if converged
                n = i
                g = conv_var(V_vec,V_before)*100
                # println("Precision of $g% reached after $n steps")
                break
            end
        end
        
        if i == N_max
            g = conv_var(V_vec,V_before)*100
            printstyled("Convergence failed after $i steps...", color=:red); println()
            println("Precision reached $g%")
        end
    end
    
    return V_vec
end

# linearised computation (for big mu)
function linearised(mu::Float64,con::Float64,
                    V0_vec::Array{Float64,1},
                    P::Array{Float64,2},
                    Nz::Int64,Nphi::Int64)
    C = I + con * 2 * mu * P
    V_vec = Symmetric(C) \ V0_vec
    V = vec_to_mat(V_vec,Nz,Nphi)
    S = dens(V,mu,con)
    
    return V, S
end

# compute the potential and surface density for a given initial potential
function compute_surf(mu::Float64,p::Float64,
                      V0_vec::Array{Float64,1},P::Array{Float64,2,},
                      con::Float64,Nz::Int64,Nphi::Int64)
    if mu <= 200.
        V_vec = conv_patch(p,V0_vec,P,mu,con)
        V = vec_to_mat(V_vec,Nz,Nphi)
        S = dens(V,mu,con)
    else
        V, S = linearised(mu,con,V0_vec,P,Nz,Nphi)
    end
    
    return V, S
end

# potential in the bulk at a polar point r for a known array of surface charge density
function V_bulk(r::Array{Float64,1},
                sigmas::Array{Float64,2},
                zs::Array{Float64,1},phis::Array{Float64,1},
                L::Float64,R::Float64)
    ab = 1.
    Nz = length(zs)
    Nphi = length(phis)
    ρ,phi,z = r
    
    V = 0.
    # sum up over all surface charges
    for i in 1:Nz
        for j in 1:Nphi
            V += sigmas[i,j] * 2π * R * L / sqrt(polar_dist(r,[R,phis[j],zs[i]],L,h)^2 + ab^2) / Nz / Nphi
        end
    end
    
    return V
end

# screening potential at an i-th dopant position
function V_screen(i::Int64,
                  pos::Array{Float64,2},
                  sigmas::Array{Float64,2},zs::Array{Float64,1},phis::Array{Float64,1},
                  L::Float64,R::Float64)
    
    # V_bulk at pos of the charge
    V = V_bulk(pos[:,i],sigmas,zs,phis,L,R)
    
    return V
end

# compute the array of screening potentials at every dopant position
function screen_pots(pos::Array{Float64,2},
                     sigmas::Array{Float64,2},zs::Array{Float64,1},phis::Array{Float64,1},
                     L::Float64,R::Float64)
    
    Vs = zeros(size(pos,2))
    
    for i in 1:size(pos,2)
        Vs[i] = V_screen(i,pos,sigmas,zs,phis,L,R)
    end
    
    return Vs
end






