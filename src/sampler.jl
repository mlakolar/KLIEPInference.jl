struct IsingSampler <: Sampleable{Multivariate,Discrete}
    p::Int64                 # number of vertices
    θ::Vector{Float64}       # parameter matrix
    nbs::Vector{Vector{Int64}}
    state::BitVector
    thin::Int64
end

Base.length(s::IsingSampler) = s.p
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
    for node=1:length(spl)
        _sample_one_coordinate!(spl, node)
    end
    spl
end

function _sample_one_coordinate!(
    spl::IsingSampler,
    node::Int64)

    θ, nbs, y = spl.θ, spl.nbs, spl.state

    v = 0.
    for nb in spl.nbs[node]
        k = trimap(node, nb)
        v += y[nb] ? θ[k] : -θ[k]
    end
    Q = exp(2. * v)
    y[node] = rand() <= Q / (Q + 1.)
end


### graph construction

# lenC --- length of a chain
# the graph consists of p / lenC chains
# sgn +1, 0, -1 --- whether edges are positive, mixed, or negative
# lb and ub are assumed positive
function chain(p::Integer, lenC::Integer=10, lb::Float64=0.1, ub::Float64=0.3, sgn::Integer=1)
    mod(p, lenC) == 0 || throw(ArgumentError("p=$p should be divisible by lenC=$lenC"))
    m = div(p*(p-1), 2)
    θ = zeros(Float64, m)

    du = Uniform(lb, ub)
    for ic=1:div(p, lenC)
        for j=1:lenC-1
            col = (ic - 1) * lenC + j
            row = col + 1

            v = rand(du)
            if sgn == -1
                v *= -1.
            elseif sgn == 0
                v *= 2. * rand(Bernoulli()) - 1.
            end
            θ[trimap(row, col)] = v
        end
    end
    θ
end


function removeEdges!(θ, numChanges::Integer=4)
    ind_change = sample(findall(!iszero, θ), numChanges; replace=false)
    for i in ind_change
        θ[i] = 0.
    end
    θ
end

function addEdges!(θ, numChanges::Integer=4, lb::Float64=0.2, ub::Float64=0.4, sgn::Integer=1)
    ind_change = sample(findall(iszero, θ), numChanges; replace=false)
    du = Uniform(lb, ub)
    for i in ind_change
        v = rand(du)
        if sgn == -1
            v *= -1.
        elseif sgn == 0
            v *= 2. * rand(Bernoulli()) - 1.
        end
        θ[i] = v
    end
    θ
end
