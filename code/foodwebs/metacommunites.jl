function sample_lonlat(N)
    lat = range(0,1,length = N)
    tuple.(rand(N), asin.(lat))
end

#constructers
"""
    metacommuntiy(coms::Array{Community})

Creates a MetaCommunity object from communities in `coms`.
"""
function metacommuntiy(coms::Array{T}) where T <: AbstractCommunity
    #get sp_vector
    sp_vec = unique(vcat([x.sp for x = coms]...))
    T_mat = [c.T for c = coms]
    loc = [c.loc for c = coms]

    #generate Sp location dictionaries
    sp_dict = Dict{UUID, Vector{Int}}()
    for i = eachindex(coms)
        for id = coms[i].ids
            if !haskey(sp_dict, id)
                sp_dict[id] = [i]
            else
                push!(sp_dict[id], i)
            end
        end
    end

    #generate distance matrix
    D = zeros(length(loc), length(loc))
    for i = eachindex(loc)
        for j = eachindex(loc)
            #calculate distance matrix
            # println(loc[i], " ", loc[j])
            if i ≠ j
                D[i,j] = acos(sin(loc[i][2]) * sin(loc[j][2]) + 
                    cos(loc[i][2]) * cos(loc[j][2]) * cos(abs(loc[i][1] - loc[j][1])))
            end
        end
    end

    #get sp vecs
    id_vec_mc = collect(keys(sp_dict))
    sp_vec_mc = hcat([sp_vec[findall([id == s.id for s = sp_vec])] for id = id_vec_mc]...)[:]

    #get communtiy Ratios
    x = [com.R for com = coms]
    @assert all( y -> y==x[1], x) "ppmr must be the same across communities"

    return MetaCommunity(coms, loc, D, T_mat, sp_vec_mc, id_vec_mc, sp_dict, coms[1].R)
end


function remove_sp_meta!(mc::MetaCommunity, a, id::UUID)
    mc.coms[a] = remove_species(mc.coms[a], id)
    
    #update loc
    filter!(x -> x != a, mc.sp_loc[id])
end

function add_sp_meta!(mc::MetaCommunity, a, sp::Species)
    mc.coms[a] = add_species(mc.coms[a], sp, calc_prod = false)
    
    #update loc
    append!(mc.sp_loc[sp.id], a)
end

"""
    move_sp_meta!(mc::MetaCommunity, a, b, id)

Moves species with `id` from community a to b in a metacommuntiy. a and b are the  1-D indexes of the communties. 
"""
function move_sp_meta!(mc::MetaCommunity, a, b, id)
    mc.coms[a], mc.coms[b] = move_species(mc.coms[a], mc.coms[b], id)

    #update loc
    append!(mc.sp_loc[id], b)
    filter!(x -> x != a, mc.sp_loc[id])
end

"""
    test_metacommunity(mc)

Tests to see if the sp_loc dictionary is correct (i.e. all species are where they should be).
"""
function check_metacommunity(mc::MetaCommunity)
    #check dicts are correct
    for (i,sp) in enumerate(mc.sp)
        coms = mc.sp_loc[sp.id]
        for c in coms
            @assert sp in mc.coms[c].sp "$sp is not in $c"
        end
    end

    for (i,c) = enumerate(mc.coms)
        for (j,id) = enumerate(c.ids)
            @assert i in mc.sp_loc[id] "Sp $id not in community $i"
        end
    end
end

"""
    get_richness(mc::MetaCommunity)

Gets vector of richness in a metacommuntiy
"""
function get_richness(mc::MetaCommunity)
    [length(c.sp) for c = mc.coms]
end