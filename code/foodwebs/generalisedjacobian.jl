"""
    structural_parameters(N::Int64, A::Matrix{Float64}, M::Vector{Float64}, n::Vector{Float64}, R::Float64, producer::Vector{Bool})

Generate structural parameters
"""
function structural_parameters(N::Int64, A::Matrix{Float64}, M::Vector{Float64}, n::Vector{Float64}, R::Float64, producer::Vector{Bool}, Tpk::Vector{Float64}, T::Float64, fp::Function = d -> exp(-abs(d)))
    #proportion of predation of j done by i
    β = A ./ (norm.(eachcol(A), 1)')
    #proportion of diet of i made up by j
    χ = A ./ (norm.(eachrow(A), 1))

    β[isnan.(β)] .= 0.0
    χ[isnan.(χ)] .= 0.0
      
    α = (M .^ -0.25) .* (fp.(T .- Tpk))


    #total growth from predation
    ρ = 1 .* .!producer
    #growth from production
    ρ̃ = 1 .- ρ
    
    #loss from predation
    σ = 1 .* (sum(A,dims=1) .!= 0)[:]
    #loss from mortality
    σ̃ = 1 .- σ

    return StructuralParameters(N,A,α,β,χ,ρ,ρ̃,σ,σ̃)
end

structural_parameters(com::Community, fp::Function = d -> exp(-abs(d))) = structural_parameters(com.N, com.A, [x.M[1] for x = com.sp], com.n, com.R, [x.producer[1] for x = com.sp], [x.Tpk for x = com.sp], com.T, fp)


function random_parameters(N::Int64, M::Int64)
    #exponent
    γ = rand(Uniform(0.5, 1.5), N, M) #[0.8, 1.5]
    λ = ones(N,N) # 1
    μ = rand(Uniform(1.0, 2.0), N, M) #[1.0, 2.0] 
    ϕ = rand(Uniform(0.0, 1.0), N, M) #[0.0, 1.0]
    ψ = rand(Uniform(0.5,1.5), N, M) #[0.5, 1.2]

    return [ExponentialParameters(γ[:,i], λ, μ[:,i], ϕ[:,i], ψ[:,i]) for i = 1:M]
end

random_parameters(com::Community, M::Int64) = random_parameters(com.N, M)

"""
    generalised_parameters(com, K)

Randomly samples a set of generalised parameters from a given interaction matrix and set of niches 
"""
function generalised_parameters(com::Community, K::Int64, fexp::Function)
    sp = structural_parameters(com)
    ep = fexp(com.N , K)

    return GeneralisedParameters.(Ref(sp), ep)
end

"""
    generalised_parameters!(com::Community, K::Int64, sp::StructuralParameters)

TBW
"""
function generalised_parameters!(com::Community, K::Int64, sp::StructuralParameters, fexp::Function)
    ep = fexp(com.N , K)
    return GeneralisedParameters.(Ref(sp), ep)
end

"""
    generalised_jacobian!(J,p::GeneralisedParameters)

Constructs the jaccobian in place as defined by the generalised model. Takes parameter values (including network structure) from the `GeneralisedParameters` object `p`.
"""
function generalised_jacobian!(J, s::StructuralParameters, e::ExponentialParameters)
    #construct generalised_jacobian
    for i = 1:s.N
        for j = 1:s.N
            if i == j #intra
                J[i,i] = s.ρ̃[i] * e.ϕ[i] + #production
                         s.ρ[i] * (e.γ[i] * s.χ[i,i] * e.λ[i,i] + e.ψ[i]) - #predation
                         s.σ̃[i] * e.μ[i] - #mortality
                         s.σ[i] * s.β[i,i] * e.ψ[i] #predation loss
                for k = 1:s.N
                    J[i,i] -= s.σ[i] * s.β[k,i] * e.λ[k,i] * ((e.γ[k] - 1) * s.χ[k,i] + 1) #predation loss
                end
            else #inter
                J[i,j] = 0

                if s.χ[i,j] != 0
                    J[i,j] = s.ρ[i] * e.γ[i] * s.χ[i,j] * e.λ[i,j] 
                end

                if s.β[j,i] != 0
                    J[i,j] -= s.σ[i] * s.β[j,i] * e.ψ[j]
                end

                for k = 1:s.N
                    if (s.β[k,i] != 0) && (s.χ[k,j] != 0)
                        J[i,j] -= s.σ[i] * (s.β[k,i] * e.λ[k,j] * (e.γ[k] - 1) * s.χ[k,j])
                    end
                end
            end

            J[i,j] *= s.α[i]
        end
    end
end

generalised_jacobian!(J,p::GeneralisedParameters) = generalised_jacobian!(J, p.s, p.e)


"""
    generalised_jacobian(p::GeneralisedParameters)

generate jacobian from parameter set p
"""
function generalised_jacobian(s::StructuralParameters, e::ExponentialParameters)
    J = zeros(s.N, s.N)
    generalised_jacobian!(J,s,e)
    return J 
end

generalised_jacobian(p::GeneralisedParameters) = generalised_jacobian(p.s, p.e)

#PSW
function proportion_stable_webs!(J::Matrix{Float64}, com::Community, fexp::Function;  N_trials::Int = 100)
    #get parameter set
    s = structural_parameters(com)
    psw = 0
    #calculate jaccobian
    p = generalised_parameters!(com, N_trials, s, fexp)
    for i = 1:N_trials
        generalised_jacobian!(J,p[i])
        λ = real(eigvals!(J)[end])
        psw += λ < 0
    end

    return psw / N_trials
end

"""
    proportion_stable_webs(com::Community, fexp::Function, N_trials::Int = 100)

Calculates the proportion stable parameter configurations for a Community.   
"""
function proportion_stable_webs(com::Community, fexp::Function; N_trials::Int = 100)
    J = zeros(com.N, com.N)
    return proportion_stable_webs!(J,com, fexp; N_trials = N_trials)
end
 
"""
    proportion_stable_webs(mc::MetaCommunity; N_trials::Int = 100)

Calculates proprtion stable webs across a metacommuntiy, returning an array of proportions. 
"""
function proportion_stable_webs(mc::MetaCommunity, fexp::Function; N_trials::Int = 100)
    props = zeros(size(mc.coms))
    for (i,c) = enumerate(mc.coms)
        if length(c.sp) > 0
            J = zeros(size(c.A))
            props[i] = proportion_stable_webs!(J,c, fexp, N_trials = N_trials)
        end
    end
    return props 
end