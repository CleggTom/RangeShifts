module FoodWebs

using Distributions
using UUIDs
using LinearAlgebra
using StatsBase
# using Graphs

include("types.jl")
include("nichemodel.jl")
include("community_utils.jl")
include("metacommunites.jl")
include("generalisedjacobian.jl")
# include("kineticmodel.jl")

end
