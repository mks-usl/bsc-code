# shift gap in order to average the chaos
function shift_gap(Es::Array{Float64,1},occ::Array{Float64,1})
    Es_occ = []
    Es_emp = []
    for i in eachindex(occ)
        if occ[i] == 1
            push!(Es_occ,Es[i])
        else
            push!(Es_emp,Es[i])
        end
    end
    middle = 0.5 * (maximum(Es_occ) + minimum(Es_emp))
    return Es .- middle
end

# average the single particle energies
function average_Es(Ndops::Int64,set::Int64,dir::String)
    
    mean_Es = zeros(Ndops*set)
    for i in 1:set
        E_adr = string("/home/usoltcev/Documents/",dir,"/energies/Es$i")
        Es = readdlm(E_adr)[:,1]
        
        O_adr = string("/home/usoltcev/Documents/",dir,"/dopantOcc/occ$i")
        Occ = readdlm(O_adr)[:,1]
        
        Es = shift_gap(Es,Occ)
        mean_Es[(i-1)*Ndops+1 : i*Ndops] += Es
    end
    return mean_Es
end

# technical function to reshape the neutral dopant position array
function shape_neut(v::Array{Any,1}) 
    l = length(v)
    npos = zeros(3,l)
    for i in eachindex(v)
        npos[:,i] = v[i]
    end
    return npos
end

# extract the positions and SPEs of the neutral dopants
function neutrals(Es_opt::Array{Float64,1},pos::Array{Float64,2},occ::Array{Float64,1},typ::Array{Float64,1})
    
    qs = @. (typ + 1) / 2 - occ
    
    neutral_acc = [] ; neutral_don = []
    Es_don = Float64[] ; Es_acc = Float64[]
    
    for i in 1:length(Es_opt)
        if qs[i] == 0
            if typ[i] == 1
                push!(neutral_don, pos[:,i])
                push!(Es_don, Es_opt[i])
            else
                push!(neutral_acc, pos[:,i])
                push!(Es_acc, Es_opt[i])
            end
        end
    end
    n_don = shape_neut(neutral_don); n_acc = shape_neut(neutral_acc)
    
    return n_don, n_acc, Es_don, Es_acc
end

# technical function to reshape the SPEs
function shape_Es(Es_neut)
    Es_arr = Float64[]
    for i in 1:length(Es_neut)
        for j in 1:length(Es_neut[i])
            push!(Es_arr,Es_neut[i][j])
        end
    end
    return Es_arr
end

# technical function to reshape the neutral dopant positions
function shape_pos(pos_neut)
    pos_arr = []
    for i in 1:length(pos_neut)
        for j in 1:size(pos_neut[i],2)
            push!(pos_arr,pos_neut[i][:,j])
        end
    end
    return shape_neut(pos_arr)
end

# function doing the disorder averaging and returning an array of averaged SPEs and the qualities of interest for neutral donors and acceptors 
function av_and_neut(Ndops::Int64,set::Int64,dir::String)
    avEs = zeros(Ndops*set)
    avEs_don = []
    avEs_acc = []
    pos_don = []
    pos_acc = []
    
    for i in 1:set # set = number of disorder configurations
        E_adr = string("/home/usoltcev/Documents/",dir,"/energies/Es$i")
        Es = readdlm(E_adr)[:,1]
        
        O_adr = string("/home/usoltcev/Documents/",dir,"/dopantOcc/occ$i")
        Occ = readdlm(O_adr)[:,1]
        
        T_adr = string("/home/usoltcev/Documents/",dir,"/dopantTyp/typ$i")
        Typ = readdlm(T_adr)[:,1]
        
        P_adr = string("/home/usoltcev/Documents/",dir,"/dopantPos/pos$i")
        Pos = readdlm(P_adr)
        
        Es = shift_gap(Es,Occ)
        avEs[(i-1)*Ndops+1 : i*Ndops] += Es
        
        n_don, n_acc, Es_don, Es_acc = neutrals(Es,Pos,Occ,Typ)
        push!(avEs_don,Es_don); push!(avEs_acc,Es_acc)
        push!(pos_don,n_don); push!(pos_acc,n_acc)
    end
    
    return avEs, shape_Es(avEs_don), shape_Es(avEs_acc), shape_pos(pos_don), shape_pos(pos_acc)
end











