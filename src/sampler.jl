struct IsingSampler <: Sampleable{Multivariate,Discrete}
    m::Int64                 # number of vertices
    θ::Matrix{Float64}       # parameter matrix
    nbs::Vector{Vector{Int64}}
    state::BitVector
    thin::Int64
end


function IsingSampler(θ::Matrix{Float64}, burn::Int64=3000, thin::Int64=1000)
    m = size(θ, 1)

    # Find neighbors
    nbs = Vector{Vector{Int64}}(undef, m)
    for node=1:m
        nbs[node] = Vector{Int64}(undef, 0)
    end
    for node=1:m
        for nb=1:m
            if node==nb
                continue
            end
            if !iszero(θ[node, nb])
                push!(nbs[node], nb)
            end
        end
    end

    # initial state
    state = rand(Bernoulli(), m) .== 1.
    spl = IsingSampler(m, θ, nbs, state, thin)
    for j=1:burn
        _sample_one!(spl)
    end
    spl
end

Base.length(s::IsingSampler) = s.m
Base.eltype(s::IsingSampler) = Bool

Distributions.rand(s::IsingSampler) = Distributions._rand!(s, BitVector(undef, length(s)))
Distributions.rand(s::IsingSampler, n::Int) = Distributions._rand!(s, BitMatrix(undef, length(s), n))

function Distributions._rand!(s::IsingSampler, x::BitMatrix)
    n = size(x, 2)
    for i=1:n
      for j=1:s.thin
        _sample_one!(s)
      end
      _sample_one!(s)
      x[:, i] .= s.state
    end
    x
end

function Distributions._rand!(s::IsingSampler, x::BitVector)
    for j=1:s.thin
        _sample_one!(s)
    end
    _sample_one!(s)
    x .= s.state
end

function _sample_one!(spl::IsingSampler)
    for node=1:spl.m
        _sample_one_coordinate!(spl, node)
    end
    spl
end

function _sample_one_coordinate!(
    spl::IsingSampler,
    node::Int64)

    θ, nbs, y = spl.θ, spl.nbs, spl.state

    v = θ[node, node]
    for nb in spl.nbs[node]
        v += y[nb] ? θ[node, nb] : -θ[node,nb]
    end
    Q = exp(2. * v)
    y[node] = rand() <= Q / (Q + 1.)
end


### graph construction

function chain(p::Int64, lb, ub)
    Θ = diagm(1 => (rand(Uniform(lb, ub), p-1) .* (2. .* rand(Bernoulli(), p-1) .- 1.)))
    Θ + Θ'
end
