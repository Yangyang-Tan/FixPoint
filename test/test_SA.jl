function du1solve(
    eqnfun,
    A0::T;
    rhomin::T=T(0.00001),
    rhomax::T=T(10),
    d::T=T(3),
    n::T=T(1),
    η::T=T(0),
    atol=1e-5,
    rtol=1e-5,
    dtmax=0.005,
    method=RadauIIA5(),
    inifun=inidu0,
) where {T}
    # Mass = T.([1.0 0 0; 0 1.0 0; 0 0 0])
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    # c1 = (2-η)*gamma(2-d/2)/(2*(4*π)^(d/2))
    p = SA[d, n, η, c1]
    f = ODEFunction(eqnfun)
    tspan = (rhomax, rhomin)
    # differential_vars = [true, true, false]
    # u0 = iniO1du1(A0, rhomax)
    u0 = inifun(A0, rhomax; d=d, n=n, η=η)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(
        prob, method; abstol=atol, reltol=rtol, maxiters=10^7, dtmax=dtmax, dtmin=1e-27
    )
    return USolution(sol)
end

#test SVector
let c1,d,η,p,f
    d=myT(2.0);η=usol_MWH[2//1, 195//100].η;n=usol_MWH[2//1, 195//100].n
    c1=(2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    p = SA[d, n, η, c1]
    f = ODEFunction(O1du1_SA)
    tspan = (usol_MWH[2//1, 195//100].rhomax, usol_MWH[2//1, 195//100].rhomin)
    u0 = SA[0.9*inidu0(usol_MWH[2//1, 195//100].A0, usol_MWH[2//1, 195//100].rhomax; d=d, n=n, η=η)...]
    prob = ODEProblem(f, u0, tspan, p)
    @time sol = solve(
        prob, RadauIIA5(); save_everystep = false,abstol=1e-23, reltol=1e-23, maxiters=10^7, dtmax=0.005, dtmin=1e-27
    )
end

#normal
let c1,d,η,p,f,sys,prob,prob_jac
    d=myT(2.0);η=usol_MWH[2//1, 195//100].η;n=usol_MWH[2//1, 195//100].n
    c1=(2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    p = [d, n, η, c1]
    f = ODEFunction(O1du1)
    tspan = (usol_MWH[2//1, 195//100].rhomax, usol_MWH[2//1, 195//100].rhomin+5.0)
    u0 = [0.99*inidu0(usol_MWH[2//1, 195//100].A0, usol_MWH[2//1, 195//100].rhomax; d=d, n=n, η=η)...]
    prob = ODEProblem(f, u0, tspan, p)
    # @time sol = solve(
    #     prob, RadauIIA5(); save_everystep = false,abstol=1e-23, reltol=1e-23, maxiters=10^7, dtmax=0.005, dtmin=1e-27
    # )
    sys=modelingtoolkitize(prob)
    prob_jac = ODEProblem(sys, [], (0.0, 1e5), jac = true)
    @time sol = solve(
        prob, TaylorMethod(150);dt=myT(0.0001),abstol=1e-23, reltol=1e-23, maxiters=10^7, dtmax=0.005, dtmin=1e-27
    )
    # prob_jac
    plot(x -> sol(x)[2], 0.1, 2.0)
end

let c1,d,η,p,f,sys,prob,prob_jac
    d=myT(2.0);η=usol_MWH[2//1, 195//100].η;n=usol_MWH[2//1, 195//100].n
    c1=(2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    p = [d, n, η, c1]
    f = ODEFunction(O1du1)
    tspan = (usol_MWH[2//1, 195//100].rhomax, usol_MWH[2//1, 195//100].rhomin)
    u0 = [0.999*inidu0(usol_MWH[2//1, 195//100].A0, usol_MWH[2//1, 195//100].rhomax; d=d, n=n, η=η)...]
    prob = ODEProblem(f, u0, tspan, p)
    # @time sol = solve(
    #     prob, RadauIIA5(); save_everystep = false,abstol=1e-23, reltol=1e-23, maxiters=10^7, dtmax=0.005, dtmin=1e-27
    # )
    sys=modelingtoolkitize(prob)
    prob_jac = ODEProblem(sys, [], (0.0, 1e5), jac = true)
    @time sol = solve(
        prob, RadauIIA5();abstol=1e-23, reltol=1e-23, maxiters=10^7, dtmax=0.005, dtmin=1e-27
    )
    # prob_jac
    plot!(x -> sol(x)[2], 0.1, 2.0)
end
using Sundials,SteadyStateDiffEq,TaylorIntegration
let diffeq,O1du1
    d=myT(2.0);η=usol_MWH[2//1, 195//100].η;n=usol_MWH[2//1, 195//100].n
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    p = [d, n, η, c1]
    @taylorize function O1du1(du, u, p, rho)
    local d =p[1]
    local n =p[2]
    local η =p[3]
    local c1 =p[4]
    du[1] = u[2]
    du[2] = u[3]
    du[3] =
        (
            (1 + u[2] + 2 * rho * u[3])^2 * ((η - 2) * u[2] + (η - 2 + d) * rho * u[3]) -
            c1 * u[3] * (3 + (n - 1) * (1 + (2 * rho * u[3]) / (1 + u[2]))^2)
        ) / (2 * c1 * rho)
    return nothing
    end
    u0 = [0.999*inidu0(usol_MWH[2//1, 195//100].A0, usol_MWH[2//1, 195//100].rhomax; d=d, n=n, η=η)...]
    tspan = (usol_MWH[2//1, 195//100].rhomax, usol_MWH[2//1, 195//100].rhomin)
    prob = ODEProblem(O1du1, u0, tspan, p)
    @time sol = solve(
        prob, DPRKN6();abstol=1e-23, reltol=1e-23, maxiters=10^7, dtmax=0.001, dtmin=1e-27
    )
    # solT = solve(prob, TaylorMethod(25), abstol=1e-21);
    # d=myT(2.0);η=usol_MWH[2//1, 195//100].η;n=usol_MWH[2//1, 195//100].n
    # c1=(2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    # p1 = [d, n, η, c1]
end
