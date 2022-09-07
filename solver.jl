function O1du1solve(
    eqnfun;
    rhomin::T,
    rho0::T,
    A0::T,
    atol=1e-5,
    rtol=1e-5,
    method=Rodas4()
) where {T}
    Mass = T.([1.0 0 0; 0 1.0 0; 0 0 0])
    f = ODEFunction(eqnfun, mass_matrix=Mass)
    tspan = (rho0, rhomin)
    # differential_vars = [true, true, false]
    u0 = iniO1du1(A0, rho0)
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, method, abstol=atol, reltol=rtol, maxiters=10^7, dtmax=0.00002, dtmin=1e-200)
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
    method=Rodas4()
) where {T}
    Mass = T.([1.0 0 0; 0 1.0 0; 0 0 0])
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    p = [d, n, η, c1]
    f = ODEFunction(eqnfun, mass_matrix=Mass)
    tspan = (rho0, rhomin)
    # differential_vars = [true, true, false]
    u0 = inidu0(A0, rho0; d=d, n=n, η=η)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, method, abstol=atol, reltol=rtol, maxiters=10^7, dtmax=0.001, dtmin=1e-400)
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
    method=Tsit5()
) where {T}
    U1 = u1fun
    U2 = u2fun
    tspan = (rho0, rhomin)
    p = (d, n, η, λ)
    u0 = iniEigen(A0, rho0, d, η, λ)
    f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, U1, U2)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, method, abstol=atol, reltol=rtol, maxiters=10^7)
    return USolution(sol)
end


###############################################
#       Optimization  functions               #
###############################################

function getA0(A0)
    sol = du0solve(
        O1du0,
        myT(A0),
        rtol=1e-20,
        atol=1e-20,
        method=RadauIIA5(),
    )
    @show A0
    return @show Nderivative(sol.t[1:100], sol.u3[1:100], sol.t[1])
end




function getA02(A0)
    sol = O1du1solve(
        O1du1,
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
    method=Tsit5()
) where {T}
    U1 = u1fun
    U2 = u2fun
    tspan = (rho0, rhomin)
    p = (d, n, η, λ)
    u0 = iniEigen(A0, rho0, d, η, λ)
    f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, U1, U2)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, method, abstol=atol, reltol=rtol, maxiters=10^7)
    return USolution(sol)
end


function getlambda0(u1fun,u2fun, rho0::T, rhomin::T, d::T=T(3), n::T=T(1), η::T=T(0))
    U1 = u1fun
    U2 = u2fun
    tspan = (rho0, rhomin)
    function lossEigen(para)
        p = (d, n, η, λ)
        f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, U1, U2)
        u0 = iniEigen(A0, rho0, d, η, λ)
        prob = ODEProblem(f, u0, tspan, p)
        sol = solve(prob, Tsit5(), p=p, saveat=tsteps)
        loss = sum(abs2, sol .- 1)
        return loss, sol
    end
    sol = solve(prob, method, abstol=atol, reltol=rtol, maxiters=10^7)
end
