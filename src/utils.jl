
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
