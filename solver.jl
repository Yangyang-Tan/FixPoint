function O1du1solve(
    eqnfun; rhomin::T, rho0::T, A0::T, atol=1e-5, rtol=1e-5, method=Rodas4()
) where {T}
    Mass = T.([1.0 0 0; 0 1.0 0; 0 0 0])
    f = ODEFunction(eqnfun; mass_matrix=Mass)
    tspan = (rho0, rhomin)
    # differential_vars = [true, true, false]
    u0 = iniO1du1(A0, rho0)
    prob = ODEProblem(f, u0, tspan)
    sol = solve(
        prob, method; abstol=atol, reltol=rtol, maxiters=10^7, dtmax=0.00002, dtmin=1e-200
    )
    return USolution(sol)
end

function du0solve(
    eqnfun,
    A0::T;
    d::T=T(3),
    n::T=T(1),
    η::T=T(0),
    rhomin::T=T(0.00001),
    rho0::T=T(10),
    atol=1e-5,
    rtol=1e-5,
    method=Rodas4(),
) where {T}
    Mass = T.([1.0 0 0; 0 1.0 0; 0 0 0])
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    p = [d, n, η, c1]
    f = ODEFunction(eqnfun; mass_matrix=Mass)
    tspan = (rho0, rhomin)
    # differential_vars = [true, true, false]
    u0 = inidu0(A0, rho0; d=d, n=n, η=η)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(
        prob, method; abstol=atol, reltol=rtol, maxiters=10^7, dtmax=0.001, dtmin=1e-400
    )
    return USolution(sol)
end

function Eigensolve(
    Eigenfun;
    u1fun,
    u2fun,
    rhomin::T,
    rho0::T,
    A0::T,
    d::T,
    n::T,
    η::T=0,
    λ::T,
    atol=1e-8,
    rtol=1e-8,
    method=Tsit5(),
) where {T}
    U1 = u1fun
    U2 = u2fun
    tspan = (rho0, rhomin)
    p = (d, n, η, λ)
    u0 = iniEigen(A0, rho0, d, η, λ)
    f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, U1, U2)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, method; abstol=atol, reltol=rtol, maxiters=10^7)
    return USolution(sol)
end

###############################################
#       Optimization  functions               #
###############################################

function getA0(A0)
    sol = du0solve(O1du0, myT(A0); rtol=1e-20, atol=1e-20, method=RadauIIA5())
    @show A0
    return @show Nderivative(sol.t[1:100], sol.u3[1:100], sol.t[1])
end

function getA02(A0)
    sol = O1du1solve(
        O1du1;
        rhomin=myT(0.00001),
        rho0=myT(10),
        A0=myT(A0),
        rtol=1e-20,
        atol=1e-20,
        method=RadauIIA5(),
    )
    @show A0
    return @show Nderivative(sol.t[1:100], sol.u3[1:100], sol.t[1])
end

function Eigensolve(
    Eigenfun;
    u1fun,
    u2fun,
    rhomin::T,
    rho0::T,
    A0::T,
    d::T,
    n::T,
    η::T=0,
    λ::T,
    atol=1e-8,
    rtol=1e-8,
    method=Tsit5(),
) where {T}
    U1 = u1fun
    U2 = u2fun
    tspan = (rho0, rhomin)
    p = (d, n, η, λ)
    u0 = iniEigen(A0, rho0, d, η, λ)
    f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, U1, U2)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, method; abstol=atol, reltol=rtol, maxiters=10^7)
    return USolution(sol)
end

function LossEigensolve(
    Eigenfun;
    u1fun,
    u2fun,
    rhomin::T,
    rho0::T,
    A0::T,
    d::T,
    n::T,
    η::T=0,
    λ::T,
    atol=1e-8,
    rtol=1e-8,
    method=Tsit5(),
) where {T}
    U1 = u1fun
    U2 = u2fun
    tspan = (rho0, rhomin)
    p = (d, n, η, λ)
    u0 = iniEigen(A0, rho0, d, η, λ)
    f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, U1, U2)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(
        prob, method; abstol=atol, reltol=rtol, maxiters=10^7, dense=false, save_on=false
    )
    return (sol.u[end] .- [1.0, dv1(λ, U1, rhomin)]) ./ [1.0, dv1(λ, U1, rhomin)]
end

function getlambda0(
    Eigenfun;
    u1fun,
    u2fun,
    rho0::T,
    rhomin::T,
    lambdalb::T,
    dlambda::T,
    d::T=T(3),
    n::T=T(1),
    η::T=T(0),
    rtol,
    atol,
    method=Tsit5(),
) where {T}
    function lossEigen(para)
        return sum(
            abs2,
            LossEigensolve(
                Eigenfun,
                u1fun=u1fun,
                u2fun=u2fun,
                rhomin=rhomin,
                rho0=rho0,
                A0=para[1],
                d=d,
                n=n,
                η=η,
                λ=para[2],
                atol=atol,
                rtol=rtol,
                method=method,
            ),
        )
    end
    # callback = function (p, l)
    #     display(l)
    #     # plt = plot(pred; ylim=(0, 6))
    #     # display(plt)
    #     # Tell Optimization.solve to not halt the optimization. If return true, then
    #     # optimization stops.
    #     return false
    # end
    u0=[-10.0,lambdalb+dlambda/2].|>T
    lb = T.([-200.0,lambdalb])
    ub = T.([200.0,lambdalb+dlambda])
    # adtype = Optimization.AutoZygote()
    # optf = Optimization.OptimizationFunction((x, p) -> lossEigen(x), Optimization.AutoForwardDiff())
    # optprob = Optimization.OptimizationProblem((x, p) -> lossEigen(x),u0, p, lb=[1.5, -200.0].|>T, ub=[1.6, -50.0].|>T)
    # result_ode = Optimization.solve(optprob, BBO_adaptive_de_rand_1_bin())
    result_ode=Optim.optimize(lossEigen,lb,ub, u0, Fminbox(NelderMead()), Optim.Options(iterations = 100,g_tol=1e-5))
    return result_ode
end

using OptimizationOptimJL
sol = solve(prob,NelderMead())
