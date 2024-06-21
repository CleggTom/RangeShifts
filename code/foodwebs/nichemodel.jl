"""
    species(C::Float64, Tpk::Vector{Float64}, N::Int64)

Method for generating multiple species at a time.
"""
function species(C::Float64, Tpk::Vector{Float64}, N::Int64; n::Vector{Float64} = rand(N))
    β = (1 / (2*C)) - 1
    r = n .* rand(Beta(1.0, β), N)
    c = rand.(Uniform.(r / 2, n))
    uuid = [UUIDs.uuid1() for i = 1:N]
    prod = [[false] for i = 1:N]
    M = [[-Inf] for i = 1:N]
    return Species.(n,r,c,Tpk,uuid,prod,M)
end

species(C::Float64, N::Int64; n::Vector{Float64} = rand(N)) = species(C, rand(N), N, n = n)


"""
    isequalcommunity(a::AbstractCommunity, b::AbstractCommunity)

Test if two community objects are equal. This forces comparison by equality (`==`) as opposed to identity ('===') which the default `isequal` method uses. 
"""
function isequalcommunity(a::Community, b::Community)
    typeof(a) != typeof(b) && return(false)
    #get keys
    com_keys = fieldnames(Community)
    #test if all fields are the same
    return getfield.(Ref(a), com_keys) == getfield.(Ref(b), com_keys) 
end

"""
    get_TL(A, N, prod, n)

calculate trophic levels from a communitiy
"""
function get_TL(com::Community)
    A, N = com.A, com.N
    χ = (A ./ (norm.(eachrow(A), 1)))
    χ[isnan.(χ)] .= 0.0
    return (inv(I(N) - χ) * ones(N , 1))[:] 
end

"""
    update_producer!(com::Community)

update the producers in the communtiy
"""
function update_producer!(com::Community)
    for (i,sp) = enumerate(com.sp)
        if (sum(com.A[i,:]) == 0) 
            sp.producer[1] = true
        end
    end    
end

"""
    check_web!(com)

Remove double and canabalistic links from a community web
"""
function check_web!(com::Community)
    A = com.A
    sp_vec = com.sp
    #remove double links
    double_links = findall(UpperTriangular(A) .* LowerTriangular(A)' .!= 0)
    #loop over and clean double links, larger consumer stays
    for link = double_links
        i,j = link.I
        if sp_vec[i].n > sp_vec[j].n
            A[j,i] = 0
        else
            A[i,j] = 0
        end
    end

    #remove canabalism
    A[diagind(A)] .= 0
end

"""
    community(sp_vec::Vector{Species}; T::Float64 = 0.5, R::Float64 = 42.0)

Generates an adjacency matrix for a given set of species using the niche model. removes spare nodes...
"""
function community(sp_vec::Vector{Species}; loc = (0.0, π / 4), T = loc[2], R::Float64 = 42.0)
    N = length(sp_vec)
    A = zeros(N,N)
    n = zeros(N)

    #loop over species
    for (i,sp_i) = enumerate(sp_vec)
        n[i] = sp_i.n
        #if producer dont add prey
        if !(sp_i.producer[1])

            for (j,sp_j) = enumerate(sp_vec)
                #if j is within range of i
                if sp_j.n < (sp_i.c + (sp_i.r/2))
                    if  sp_j.n > (sp_i.c - (sp_i.r/2))
                            A[i,j] = 1
                    end
                end
            end

        end
    end

    #ids
    ids = [x.id for x = sp_vec]

    #create inital community with TL = 0
    com = Community(N, A, sp_vec, ids, zeros(N), n, T, loc, R)

    #remove double and canabalistic links
    check_web!(com)

    #remove isolated nodes
    to_remove = ids[findall((sum(com.A,dims = 1)[:] .== 0) .& (sum(com.A,dims = 2)[:] .== 0))]

    for rm_indx = to_remove
        com = remove_species(com, rm_indx)
    end
    
    #calculate producer
    update_producer!(com)
    # #calculate TL
    com.TL .= get_TL(com)
    #calculate Mass
    [sp.M[1] = com.R ^ com.TL[i] for (i,sp) = enumerate(com.sp)]

    return com
end