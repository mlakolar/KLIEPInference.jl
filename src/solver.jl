abstract type KLIEPSolver end

struct CD_KLIEP <: KLIEPSolver end
struct SCS_KLIEP <: KLIEPSolver end
struct Mosek_KLIEP <: KLIEPSolver end


KLIEP(Ψx, Ψy, ::SCS_KLIEP) = _KLIEP(Ψx, Ψy, JuMP.with_optimizer(SCS.Optimizer, verbose=1))
KLIEP(Ψx, Ψy, ::Mosek_KLIEP) = _KLIEP(Ψx, Ψy, JuMP.with_optimizer(MosekOptimizer, QUIET=false)

function _KLIEP(Ψx, Ψy, solver)
    m, nx = size(Ψx)
    ny = size(Ψy, 2)
    @assert m == size(Ψy, 1)

    lny = log(ny)

    problem = JuMP.Model(solver)
    JuMP.@variable(problem, Θ[1:m])
    JuMP.@variable(problem, t)
    JuMP.@variable(problem, u[1:ny])

    JuMP.@constraint(problem, sum(u) <= 1.)
    for j=1:ny
        JuMP.@constraint(problem, [dot(Θ, Ψy[:, j]) - t - lny, 1., u[j]] in MOI.ExponentialCone())
    end

    JuMP.@objective(problem, Min, -sum(Ψx'*Θ) / nx + t)
    JuMP.optimize!(problem)

    JuMP.result_value.(Θ)
end
