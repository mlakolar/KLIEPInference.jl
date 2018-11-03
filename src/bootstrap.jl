struct BootstrapEstimates
    θhat::Vector{Float64}
    θb::Matrix{Float64}     # each column represents a bootstrap estimate
end


function boot_KLIEP(Ψx, Ψy; bootSamples::Int64=300)

    θ = KLIEP(Ψx, Ψy, CD_KLIEP())
    θhat = convert(Vector, θ)
    θb = Matrix{Float64}(undef, length(θhat), bootSamples)

    m, nx = size(Ψx)
    ny = size(Ψy, 2)
    bΨx = similar(Ψx)
    bΨy = similar(Ψy)
    x_ind = Vector{Int64}(undef, nx)
    y_ind = Vector{Int64}(undef, ny)

    for b=1:bootSamples
       sample!(1:nx, x_ind)
       sample!(1:ny, y_ind)

       _fill_boot_Psi!(bΨx, Ψx, x_ind)
       _fill_boot_Psi!(bΨy, Ψy, y_ind)

       KLIEP!(θ, bΨx, bΨy)
       θb[:, b] .= θ
    end

    BootstrapEstimates(θhat, θb)
end

function simulCI(straps::BootstrapEstimates, α::Float64=0.95)
    m, bootSamples = size(straps.θb)

    infNormDist = Vector{Float64}(undef, bootSamples)
    CI = Matrix{Float64}(undef, m, 2)

    for b=1:bootSamples
        infNormDist[b] = norm_diff(straps.θhat, view(straps.θb, :, b), Inf)
    end
    x = quantile!(infNormDist, 0.95)
    CI[:, 1] .= straps.θhat .- x
    CI[:, 2] .= straps.θhat .+ x

    CI
end

function simulCIstudentized(straps::BootstrapEstimates, α::Float64=0.95)
    m, bootSamples = size(straps.θb)

    infNormDist = Vector{Float64}(undef, bootSamples)
    w = reshape(std(straps.θb; dims = 2, corrected = false), :)

    CI = Matrix{Float64}(undef, m, 2)
    tmp = Vector{Float64}(undef, m)

    for b=1:bootSamples
        tmp .= (straps.θhat .- straps.θb[:, b]) ./ w
        infNormDist[b] = norm(tmp, Inf)
    end
    x = quantile!(infNormDist, 0.95)
    @. CI[:, 1] = straps.θhat - x * w
    @. CI[:, 2] = straps.θhat + x * w

    CI
end


function _fill_boot_Psi!(bΨ, Ψ, b_ind)
  m, n = size(Ψ)
  for col=1:n
      for row=1:m
          bΨ[row, col] = Ψ[row, b_ind[col]]
      end
  end
end
