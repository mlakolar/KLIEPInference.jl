
function Ψising(X)
    m, n = size(X)

    p = div(m * (m - 1), 2)
    out = zeros(Float64, p, n)

    for i=1:n
        ind = 0
        for col=1:m
            for row=col+1:m
                ind += 1
                out[ind, i] = xor(X[col, i], X[row, i]) ? -1. : 1.
            end
        end
    end
    out
end

function KLIEP_Hessian(θ, Ψy)
  m, n = size(Ψy)

  w = transpose(Ψy) * θ
  w .= exp.(w)
  w ./= mean(w)

  StatsBase.cov(Ψy, weights(w), 2; corrected=false)
end



function _maxabs(a, b)
    maximum(x -> abs(x[1]-x[2]), zip(a, b))
end
