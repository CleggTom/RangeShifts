print("nthreads: ", Threads.nthreads())

using Pkg

Pkg.activate(".")

using Distributions
using UUIDs
using LinearAlgebra
using StatsBase
using CairoMakie
using JLD2

include("./foodwebs/FoodWebs.jl")
fw = FoodWebs

#creating a stable metacommunity
function get_stable_com(sp_vec, lonlat, θ = 0.75, fe::Function = fw.random_parameters)
    psw = 0.0
    com = 0
    while psw < θ
        sp_filter = sp_vec[fw.sample_community(30, lonlat[2], sp_vec, d -> exp(-10abs(d)))]
        com = fw.community(sp_filter, loc = lonlat)
        if com.N > 20
            psw = fw.proportion_stable_webs(com, fw.random_parameters, N_trials = 100)
        end
    end
    return(com)
end

function get_stable_metacom(sp_vec, lonlat, θ = 0.7, fe::Function = fw.random_parameters)
    coms = Vector{fw.Community}(undef, length(lonlat))
    for t = 1:length(lonlat)
        coms[t] = get_stable_com(sp_vec, lonlat[t], θ, fe)
    end
    
    return(fw.metacommuntiy(coms))
end

function get_disperal_lists(mc, ΔT; fpd::Function = d -> exp(-10abs(d)), αd = 0.75)
    #find maximum distance
    M = vcat([s.M for s = mc.sp]...)
    d = M .^ αd
    dmax =  d ./ maximum(d)
    
    #loop through species and find dispersal
    from = Vector{Vector{Int}}(undef, length(mc.sp))
    to = Vector{Vector{Int}}(undef, length(mc.sp))
    
    for i = 1:length(mc.sp)
        sp = mc.sp[i]

        if !sp.producer[1] 
            
            #caluculate distances
            reachable = dmax[i] .> mc.D[mc.sp_loc[sp.id],:]
            
            pd_after = fpd.(mc.T_mat .- sp.Tpk .- ΔT) 
            
            #combine to get probabilities to test
            to_test = vcat([findall(r) for r = eachrow(reachable)]...)
            to_disperse = rand(length(to_test)) .< pd_after[to_test]
            from[i] = deepcopy(mc.sp_loc[sp.id])
            to[i] = unique(to_test[to_disperse])
        else
            from[i] = deepcopy(mc.sp_loc[sp.id])
            to[i] = []
        end
    end

    return(from,to)
end

function move_and_remove!(mc, from, to)
    # #move and remove species
    for i = 1:length(mc.sp)
        #add species  
        for t = to[i]
            if t ∉ from[i]
                fw.add_sp_meta!(mc, t, mc.sp[i])
            end
        end
        #remove species
        for f = from[i]
            if f ∉ to[i]
                fw.remove_sp_meta!(mc, f, mc.sp[i].id)
            end
        end
    end
end

function exp_parameters(N::Int64, M::Int64)
    #exponent
    γ = rand(Uniform(0.8, 1.5), N, M) #[0.8, 1.5]
    λ = ones(N,N) # 1
    μ = rand(Uniform(1.0, 2.0), N, M) #[1.0, 2.0] 
    ϕ = rand(Uniform(0.0, 1.0), N, M) #[0.0, 1.0]
    ψ = rand(Uniform(0.5,1.2), N, M) #[0.5, 1.2]

    return [fw.ExponentialParameters(γ[:,i], λ, μ[:,i], ϕ[:,i], ψ[:,i]) for i = 1:M]
end

function mean_M(com)
    if com.N > 0
        return mean([s.M[1] for s = com.sp])
    else
        return 0
    end
end

Nrep = 500
NT = 50
Nα = range(-0.75,0.75, length = 5)
NΔT = range(0.0, 0.5, length = 5)

results_psw = zeros(Nrep, 2, 5, 5, NT)
results_N = zeros(Nrep, 2, 5, 5, NT)
results_M = zeros(Nrep, 2, 5, 5, NT)

# itt = [0]
for (i,v) = enumerate(Nα)
    for (j,w) = enumerate(NΔT)
        itt = [0]
        Npool = 500
        Tpk = rand(Npool) * (π/2)
        sp_vec = fw.species(0.1,Tpk,Npool)
        
        Threads.@threads for k = 1:Nrep
                itt[1] += 1
                # print("\r",itt[1])
                print("\r $i  $j ", itt[1])
                lonlat = fw.sample_lonlat(NT)

                mc = get_stable_metacom(sp_vec, lonlat, 0.75)
            
                results_psw[k,1,i,j,:] .= fw.proportion_stable_webs(mc, exp_parameters, N_trials = 100)
                results_N[k,1,i,j,:] .= [c.N for c = mc.coms]
                results_M[k,1,i,j,:] .= mean_M.(mc.coms)
        
                f,t = get_disperal_lists(mc, w, αd = v)
                try
                    move_and_remove!(mc, f, t)
                    results_psw[k,2,i,j,:] .= fw.proportion_stable_webs(mc, exp_parameters, N_trials = 100)
                    results_N[k,2,i,j,:] .= [c.N for c = mc.coms]
                    results_M[k,2,i,j,:] .= mean_M.(mc.coms)
                catch
                    results_psw[k,2,i,j,:] .= 0
                    results_N[k,2,i,j,:] .= 0
                    results_M[k,2,i,j,:] .= 0
                end

        end
    end
end

save("./results/results.jld2", Dict("psw" => results_psw, "N" => results_N, "M" => results_M))