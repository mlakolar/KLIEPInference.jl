
function rhat(θ, Ψy)
    r = transpose(Ψy) * θ
    r .= exp.(r)
    r ./ mean(r)
end

KLIEP_Hessian(θ, Ψy) = StatsBase.cov(Ψy, weights(rhat(θ, Ψy)), 2; corrected=false)

function Ψising(X)
    m, n = size(X)
    p = div(m * (m - 1), 2)
    out = zeros(Float64, p, n)
    for i = 1:n
        ind = 0
        for row = 2:m
            for col = 1:row-1
                ind += 1
                out[ind, i] = xor(X[col, i], X[row, i]) ? -1. : 1.
            end
        end
    end
    out
end

function _find_supp3(θ, ω, k)
    supp1 = findall(!iszero, θ)
    supp2 = findall(!iszero, ω)
    supp3 = union(supp1, supp2)
    pos_k = findfirst(isequal(k), supp3)
    if pos_k == nothing
        push!(supp3, k)
    else
        supp3[pos_k], supp3[end] = supp3[end], supp3[pos_k]
    end
    supp3
end

function _maxabs(a, b)
    maximum(x -> abs(x[1]-x[2]), zip(a, b))
end

# maps (i, j) coordinate of a symmetric matrix to a packed format
# diagonal of the matrix is not included
# lower triangular part is mapped row by row
function trimap(i::Integer, j::Integer)
    if i < j
        trimap(j, i)
    else
        j + div((i - 1) * (i - 2), 2)
    end
end

# mapping from packed storage to (i, j)
function itrimap(k::Integer)
    i = convert(Int, ceil((1. + sqrt(1. + 8. * k)) / 2.))
    j = k - div((i - 1) * (i - 2), 2)
    CartesianIndex(i, j)
end

# vector -> symmetric matrix
function unpack(θ)
    m = length(θ)
    p = convert(Int, ceil((1. + sqrt(1. + 8. * m)) / 2.))
    out = zeros(Float64, p, p)
    for i=1:m
        out[itrimap(i)] = θ[i]
    end
    out + out'
end

# symmetric matrix -> vector
function pack(Θ)
    m = size(Θ, 1)
    p = div(m * (m - 1), 2)
    out = zeros(Float64, p)
    for j = 1:(m-1)
        for i = (j+1):m
            out[trimap(i,j)] = Θ[i,j]
        end
    end
    out
end


####################################
#
# loss En[-θ'Ψx(i)] + log( En[exp(θ'Ψy(i))] )
#
####################################
struct CDKLIEPLoss <: CoordinateDifferentiableFunction
    μx::Vector{Float64}   # mean(Ψx, dims = 2)
    Ψy::Matrix{Float64}
    r::Vector{Float64}    # ny dim vector --- stores θ'Ψy(i)
    p::Int64
end

function CDKLIEPLoss(Ψx::Matrix{Float64}, Ψy::Matrix{Float64})
    (p = size(Ψx, 1)) == size(Ψy, 1) || throw(DimensionMismatch())

    CDKLIEPLoss(vec(mean(Ψx, dims=2)), Ψy, zeros(size(Ψy, 2)), p)
end

CoordinateDescent.numCoordinates(f::CDKLIEPLoss) = f.p

function CoordinateDescent.initialize!(f::CDKLIEPLoss, x::SparseIterate)

    mul!(f.r, transpose(f.Ψy), x)

    nothing
end

function CoordinateDescent.gradient(
    f::CDKLIEPLoss,
    x::SparseIterate,
    j::Int64)

    μx, Ψy, r = f.μx, f.Ψy, f.r

    -μx[j] + mean( exp.(r) .* Ψy[j, :] ) / mean(exp, r)
end

function CoordinateDescent.descendCoordinate!(
    f::CDKLIEPLoss,
    g::Union{ProxL1, ProxZero},
    x::SparseIterate,
    j::Int64)

    μx, Ψy, r = f.μx, f.Ψy, f.r

    mean_exp_r = mean(exp, r)
    mean_exp_r_Ψy = mean( exp.(r) .* Ψy[j, :] )
    mean_exp_r_Ψy2 = mean( exp.(r) .* Ψy[j, :].^2. )

    @inbounds grad = -μx[j] + mean_exp_r_Ψy / mean_exp_r
    H = mean_exp_r_Ψy2 / mean_exp_r - mean_exp_r_Ψy^2 / mean_exp_r^2.

    @inbounds oldVal = x[j]
    @inbounds x[j] -= grad / H
    newVal = cdprox!(g, x, j, 1. / H)
    h = newVal - oldVal

    # update internals
    for i=1:length(r)
        @inbounds r[i] += h*Ψy[j, i]
    end
    h
end
