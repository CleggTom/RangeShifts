#Community generation
"""
    sample_community(N::Int64, T::Float64, sp_vec::Vector{Species}, fp::Function = d -> exp(-20abs(d)))

Sample a community of a given size and temperature from a global species list. Probabiltiy of selecting species is based on function fp. 
"""
function sample_community(N::Int64, T::Float64, sp_vec::Vector{Species}, fp::Function = d -> exp(-20abs(d)))
    sample(1:length(sp_vec), Weights(fp.([sp.Tpk for sp = sp_vec] .- T)), N, replace = false)
end


"""
    add_species(com::Community, sp::Species)

Add species `sp` to a community and partially construct the niche web. Will assume species is not a producer
"""
function add_species(com::Community, sp::Species; calc_prod = true)
    @assert !in(sp.id, com.ids) "sp already present"
    # get new sp list
    sp_vec_new = vcat(com.sp, [sp])
    n_new = vcat(com.n, [sp.n])

    #ids
    ids = [x.id for x = sp_vec_new]

    #add sp without recalculating whole web structure
    
    #new array
    A = zeros(size(com.A) .+ 1) 
    A[1:end-1 , 1:end-1] .= com.A
    #loop over sp
    for i = 1:size(com.A)[1] 
        #consumption
        if !sp.producer[1]
            if (com.sp[i].n < sp.c + (sp.r/2) && com.sp[i].n > sp.c - (sp.r/2))
                A[end,i] = 1
            end
        end

        #predation
        if !com.sp[i].producer[1]
            if (sp.n < com.sp[i].c + (com.sp[i].r/2) && sp.n > com.sp[i].c - (com.sp[i].r/2))
                A[i,end] = 1
            end
        end
    end

    if calc_prod
        if sum(A[end,:]) == 0
            sp_vec_new[end].producer[1] = true
        end
    end

    TL = get_TL(com)

    new_com = Community(com.N + 1 , A, sp_vec_new, ids, TL,  n_new, com.T, com.loc, com.R)
    
    check_web!(new_com)

    return new_com
end


"""
    remove_species(com::Community, id::UUID)

Remove species identified with `id`.
"""
function remove_species(com::Community, id::UUID)
    @assert id in com.ids "id not in community"
    #get all indicies that dont match (to keep)
    indx = com.ids .!= id
    #get new species, id and niche value lists
    n_sp = com.sp[indx]
    n_ids = com.ids[indx]
    n_n = com.n[indx]
    #create new community 
    new_com = Community(com.N - 1, com.A[indx,indx], n_sp, n_ids , zeros(com.N - 1), n_n, com.T, com.loc, com.R)

    return new_com
end

# """
#     remove_species(com::Community, sp::Species)

# Remove species `sp`
# """
# function remove_species(com::Community, sp::Species)
#     return remove_species(com, sp.id)
# end

"""
    move_species(com1::Community,com2::Community,id::UUID)

Move species identified by `id` from `com1` to `com2`
"""
function move_species(com1::Community,com2::Community, id::UUID)
    @assert id in com1.ids "no species with id in com1"
    @assert !(id in com2.ids) "id aleady in com2"
    #add to community 2
    com2_new = add_species(com2, com1.sp[findfirst(com1.ids .== id)], calc_prod = false)

    #remove from community 1
    com1_new = remove_species(com1, id)

    return com1_new, com2_new
end
