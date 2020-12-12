function tuneλ(λ::Vector{Float64}, Ψx, Ψy, idx::Int64, a::Float64=0.9, nrep::Int64=1000)
    nλ = length(λ)

    nx = size(Ψx, 2)
    p, ny = size(Ψy)

    θ0 = zeros(Float64, 2, nλ)
    for t = 1:nλ
        θ = spKLIEP(Ψx, Ψy, λ[t], CD_KLIEP())
        supp = sort(union(idx, findall(!iszero, θ)))
        spKLIEP_refit!(θ, Ψx, Ψy, supp)

        H = KLIEP_Hessian(θ, Ψy)
        ω = Hinv_row(H, idx, round(sqrt(2. * log(p) / ny), digits=3))
        supp = sort(union(idx, findall(!iszero, ω)))
        Hinv_row_refit!(ω, H, idx, supp)

        θ0[:,t] = debias_KLIEP(idx, θ, ω, Ψx, Ψy)
    end

    mx = round(Int,a*nx)
    my = round(Int,a*ny)

    λ1 = sqrt(1/a) * λ
    λ2 = round(sqrt(2. * log(p) / my), digits=3)

    θhat = zeros(Float64, 2, nλ, nrep)
    for rep = 1:nrep
        if mod(rep, div(nrep, 10))==0
            @printf("%d / %d\n", rep, nrep)
        end

        Ψxsub = Ψx[:,sample(1:nx, mx, replace=false)]
        Ψysub = Ψy[:,sample(1:ny, my, replace=false)]

        for t = 1:nλ
            θ = spKLIEP(Ψxsub, Ψysub, λ1[t], CD_KLIEP())
            supp = sort(union(idx, findall(!iszero, θ)))
            spKLIEP_refit!(θ, Ψxsub, Ψysub, supp)

            H = KLIEP_Hessian(θ, Ψysub)
            ω = Hinv_row(H, idx, λ2)
            supp = sort(union(idx, findall(!iszero, ω)))
            Hinv_row_refit!(ω, H, idx, supp)

            θhat[:,t,rep] = KLIEP_debias(idx, θ, ω, Ψxsub, Ψysub)
        end
    end

    θhat[1,:,:] = broadcast(/, broadcast(-, θhat[1,:,:], θ0[1,:]), std(θhat[1,:,:],dims=2))
    θhat[2,:,:] = broadcast(/, broadcast(-, θhat[2,:,:], θ0[2,:]), std(θhat[2,:,:],dims=2))

    KSstats = zeros(Float64, 2, nλ)
    for t = 1:nλ
        KSstats[1,t] = HypothesisTests.ksstats(θhat[1,t,:], Normal())[2]
        KSstats[2,t] = HypothesisTests.ksstats(θhat[2,t,:], Normal())[2]
    end

    KSstats
end
