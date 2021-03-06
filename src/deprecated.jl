
function spKLIEP_refit!(x::SparseIterate, Ψx, Ψy, supp::Vector{Int64})
    w = ones(Float64, length(x)) * 1e10

    for k in supp
        w[k] = 0.
    end

    f = CDKLIEPLoss(Ψx, Ψy)
    g = ProxL1(1., w)
    coordinateDescent!(x, f, g)
end

spKLIEP_refit!(x::SparseIterate, Ψx, Ψy) = spKLIEP_refit!(x::SparseIterate, Ψx, Ψy, findall(!iszero, x))



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

function boot_spKLIEP(Ψx, Ψy, θfs, Hinv; bootSamples::Int64=300)

    # generate bootstrap samples first since we want to use the same
    # samples for each coordinate j=1...m
    m, nx = size(Ψx)
    ny = size(Ψy, 2)

    x_ind = Matrix{Int16}(undef, nx, bootSamples)
    y_ind = Matrix{Int16}(undef, ny, bootSamples)

    for b=1:bootSamples
       sample!(1:nx, view(x_ind, :, b))
       sample!(1:ny, view(y_ind, :, b))
    end
    θhat = Vector{Float64}(undef, m)
    θb = Matrix{Float64}(undef, length(θhat), bootSamples)

    for j=1:m
        supp3 = findSupp3(θfs, Hinv[j], j)

        # obtain θhat_j
        bΨx = Ψx[supp3, :]
        bΨy = Ψy[supp3, :]
        bbΨx = similar(bΨx)
        bbΨy = similar(bΨy)
        θ = KLIEP(bΨx, bΨy, CD_KLIEP())
        θhat[j] = θ[end]

        # bootstrap
        for b=1:bootSamples
           _fill_boot_Psi!(bbΨx, bΨx, view(x_ind, :, b))
           _fill_boot_Psi!(bbΨy, bΨy, view(y_ind, :, b))

           KLIEP!(θ, bbΨx, bbΨy)
           θb[j, b] = θ[end]
        end
    end

    BootstrapEstimates(θhat, θb)
end

function boot_spKLIEPfull(Ψx, Ψy, θfs, Hinv, λ1, λ2; bootSamples::Int64=300)

    m, nx = size(Ψx)
    ny = size(Ψy, 2)

    θhat = Vector{Float64}(undef, m)
    θb = Matrix{Float64}(undef, length(θhat), bootSamples)

    # compute third stage estimator first
    for j=1:m
        supp3 = findSupp3(θfs, Hinv[j], j)

        # obtain θhat_j
        bΨx = Ψx[supp3, :]
        bΨy = Ψy[supp3, :]
        θ = KLIEP(bΨx, bΨy, CD_KLIEP())
        θhat[j] = θ[end]
    end

    x_ind = Vector{Int64}(undef, nx)
    y_ind = Vector{Int64}(undef, ny)
    Ψx_boot = similar(Ψx)
    Ψy_boot = similar(Ψy)

    for b=1:bootSamples
        sample!(1:nx, x_ind)
        sample!(1:ny, y_ind)

        _fill_boot_Psi!(Ψx_boot, Ψx, x_ind)
        _fill_boot_Psi!(Ψy_boot, Ψy, y_ind)

        # first step
        θfs_boot = deepcopy(θfs)
        spKLIEP!(θfs_boot, Ψx_boot, Ψy_boot, λ1, CD_KLIEP())
        spKLIEP_refit!(θfs_boot, Ψx_boot, Ψy_boot)

        # second and third step
        H_boot = KLIEP_Hessian(θfs_boot, Ψy_boot)
        for j=1:m
            ω = Hinv_row(H_boot, j, λ2)
            supp3_boot = findSupp3(θfs_boot, ω, j)

            bΨx_boot = Ψx_boot[supp3_boot, :]
            bΨy_boot = Ψy_boot[supp3_boot, :]

            θ = KLIEP(bΨx_boot, bΨy_boot, CD_KLIEP())
            θb[j, b] = θ[end]
        end
    end

    BootstrapEstimates(θhat, θb)
end

# S_delta --- support of Delta
function boot_oracleKLIEP(Ψx, Ψy, S_delta; bootSamples::Int64=300)

    # generate bootstrap samples first since we want to use the same
    # samples for each coordinate j=1...m
    m, nx = size(Ψx)
    ny = size(Ψy, 2)

    x_ind = Matrix{Int16}(undef, nx, bootSamples)
    y_ind = Matrix{Int16}(undef, ny, bootSamples)

    for b=1:bootSamples
       sample!(1:nx, view(x_ind, :, b))
       sample!(1:ny, view(y_ind, :, b))
    end
    θhat = Vector{Float64}(undef, m)
    θb = Matrix{Float64}(undef, length(θhat), bootSamples)

    for j=1:m
        S = copy(S_delta)
        pos_j = findfirst(isequal(j), S)
        if pos_j == nothing
            push!(S, j)
        else
            S[pos_j], S[end] = S[end], S[pos_j]
        end

        # obtain θhat_j
        bΨx = Ψx[S, :]
        bΨy = Ψy[S, :]
        bbΨx = similar(bΨx)
        bbΨy = similar(bΨy)
        θ = KLIEP(bΨx, bΨy, CD_KLIEP())
        θhat[j] = θ[end]

        # bootstrap
        for b=1:bootSamples
           _fill_boot_Psi!(bbΨx, bΨx, view(x_ind, :, b))
           _fill_boot_Psi!(bbΨy, bΨy, view(y_ind, :, b))

           KLIEP!(θ, bbΨx, bbΨy)
           θb[j, b] = θ[end]
        end
    end

    BootstrapEstimates(θhat, θb)
end

function boot_gaussKLIEP(Ψx, Ψy, θhat, Hinv; bootSamples::Int64=300)

    # generate bootstrap samples first since we want to use the same
    # samples for each coordinate j=1...m
    m, nx = size(Ψx)
    ny = size(Ψy, 2)

    # compute weights corresponding to rhat
    w = transpose(Ψy) * θhat
    w .= exp.(w)
    w ./= mean(w)

    # means
    mux = vec( mean(Ψx, dims = 2) )
    muy = vec( mean(Ψy, weights(w), 2) )

    θb = Matrix{Float64}(undef, length(θhat), bootSamples)

    # store Gaussian multipliers
    gx = Vector{Float64}(undef, nx)
    gy = Vector{Float64}(undef, ny)

    tmp1 = Vector{Float64}(undef, m)
    tmp2 = Vector{Float64}(undef, m)
    # bootstrap
    for b=1:bootSamples
        rand!(Normal(), gx)
        rand!(Normal(), gy)

        fill!(tmp1, 0.)
        fill!(tmp2, 0.)
        for i=1:nx
            @. tmp1 += (Ψx[:, i] - mux) * gx[i]
        end

        for i=1:ny
            @. tmp2 += (Ψy[:, i] * w[i] - muy) * gy[i]
        end
        @. tmp1 = tmp1 / nx + tmp2 / ny

        for j=1:m
            θb[j, b] = θhat[j] + dot(Hinv[j], tmp1)
        end
    end

    BootstrapEstimates(θhat, θb)
end
