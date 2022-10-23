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
    method=Rodas4(),
) where {T}
    # Mass = T.([1.0 0 0; 0 1.0 0; 0 0 0])
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    p = [d, n, η, c1]
    f = ODEFunction(eqnfun)
    tspan = (rhomax, rhomin)
    # differential_vars = [true, true, false]
    # u0 = iniO1du1(A0, rhomax)
    u0 = inidu0(A0, rhomax; d=d, n=n, η=η)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(
        prob, method; abstol=atol, reltol=rtol, maxiters=10^7, dtmax=dtmax, dtmin=1e-27
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
    rhomax::T=T(10),
    atol=1e-5,
    rtol=1e-5,
    method=Rodas4(),
) where {T}
    Mass = T.([1.0 0 0; 0 1.0 0; 0 0 0])
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    p = [d, n, η, c1]
    f = ODEFunction(eqnfun; mass_matrix=Mass)
    tspan = (rhomax, rhomin)
    # differential_vars = [true, true, false]
    u0 = inidu0(A0, rhomax; d=d, n=n, η=η)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(
        prob, method; abstol=atol, reltol=rtol, maxiters=10^7, dtmax=0.0001, dtmin=1e-400
    )
    return USolution(sol)
end
function Eigensolve(
    Eigenfun;
    u1fun,
    u2fun,
    rhomin::T,
    rhomax::T,
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
    tspan = (rhomax, rhomin)
    p = (d, n, η, λ)
    u0 = iniEigen(A0, rhomax, d, η, λ)
    f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, U1, U2)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, method; abstol=atol, reltol=rtol, maxiters=10^7, dtmax=0.005)
    return USolution(sol)
end

###############################################
#       Optimization  functions               #
###############################################

# function getA0(A0)
#     sol = du0solve(O1du0, myT(A0); rtol=1e-20, atol=1e-20, method=RadauIIA5())
#     @show A0
#     return @show Nderivative(sol.t[1:100], sol.u3[1:100], sol.t[1])
# end

function getA0(A0::T; d::T=T(3.0), n::T=T(1.0), η::T=T(0.0)) where {T}
    sol = du0solve(
        O1du0, myT(A0); rtol=1e-20, atol=1e-20, method=RadauIIA5(), d=d, n=n, η=η
    )
    @show A0
    return @show Nderivative(sol.t[1:100], sol.u3[1:100], sol.t[1])
end

function getA02(
    A0::T;
    rhomin::T=T(10^-8),
    rhomax::T=T(10),
    d::T=T(3.0),
    n::T=T(1.0),
    η::T=T(0.0),
    rtol=10^-22,
    atol=10^-22,
    dtmax=0.005,
) where {T}
    sol2 = du1solve(
        O1du1,
        A0;
        rhomin=rhomin,
        rhomax=rhomax,
        rtol=rtol,
        atol=atol,
        method=RadauIIA5(),
        d=d,
        n=n,
        η=η,
        dtmax=dtmax,
    )
    println(A0)
    # return @show Nderivative(sol2.t[1:100], sol2.u3[1:100], sol2.t[1])
    return @show sol2.sol(sol2.t[1], Val{1})[3]
end

function LossEigensolve(
    Eigenfun;
    u1fun,
    u2fun,
    rhomin::T,
    rhomax::T,
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
    tspan = (rhomax, rhomin)
    p = (d, n, η, λ)
    u0 = iniEigen(A0, rhomax, d, η, λ)
    f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, U1, U2)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(
        prob, method; abstol=atol, reltol=rtol, maxiters=10^7, dense=false, save_on=false
    )
    # return (sol.u[end] .- [1.0, dv1(λ, U1, rhomin)]) ./ [1.0, dv1(λ, U1, rhomin)]
    v0 = sol.u[end][1]
    v1 = sol.u[end][2]
    v1real = dv1(; λ=λ, v0=v0, U1=U1, rhomin=rhomin, d=d, n=n, η=η)
    return (v1 - v1real)^2
end

function LossEigensolve2(
    Eigenfun;
    u1fun,
    u2fun,
    rhomin::T,
    rhomax::T,
    A0::T,
    d::T,
    n::T,
    η::T,
    λ::T,
    atol=1e-8,
    rtol=1e-8,
    method=Tsit5(),
) where {T}
    U1 = u1fun
    U2 = u2fun
    tspan = (rhomax, rhomin)
    p = (d, n, η, λ)
    u0 = iniEigen(A0, rhomax, d, η, λ)
    f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, U1, U2)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(
        prob, method; abstol=atol, reltol=rtol, maxiters=10^7, dtmax=0.005,dense=false, save_on=false
    )
    # return (sol.u[end] .- [1.0, dv1(λ, U1, rhomin)]) ./ [1.0, dv1(λ, U1, rhomin)]
    v0 = sol.u[end][1]
    v1 = sol.u[end][2]
    v1real = dv1(; λ=λ, v0=v0, U1=U1, rhomin=rhomin, d=d, n=n, η=η)
    return (v1 - v1real)
end

function findUmin(;
    lb::T,
    ub::T,
    rhomin::T,
    rhomax::T,
    d::T,
    n::T,
    η::T,
    tol_usol=1e-8,
    tol_A0=1e-15,
    method=Brent(),
    dtmax=0.005,
) where {T}
    function fun(x)
        return getA02(
            x; rhomin=rhomin, rhomax=rhomax, d=d, n=n, η=η, rtol=tol_usol, atol=tol_usol,dtmax=dtmax
        )^2
    end
    opt1 = Optim.optimize(fun, lb, ub, method; rel_tol=T(tol_A0), abs_tol=T(tol_A0))
    @show optmin1 = opt1.minimizer
    sol = du1solve(
        O1du1,
        optmin1;
        rhomin=rhomin,
        rhomax=rhomax,
        rtol=tol_usol,
        atol=tol_usol,
        d=d,
        n=n,
        η=η,
        dtmax=dtmax,
        method=RadauIIA5(),
    )
    return USolution(sol.sol; A0=optmin1, rhomin=rhomin, rhomax=rhomax, d=d, n=n, η=η)
end

function findUmin2(;
    lb::T,
    ub::T,
    rhomin::T,
    rhomax::T,
    d::T,
    n::T,
    η::T,
    tol_usol=1e-8,
    tol_A0=1e-15,
    method=Brent(),
    dtmax=0.005,
) where {T}
    function fun(x)
        return getA02(
            x; rhomin=rhomin, rhomax=rhomax, d=d, n=n, η=η, rtol=tol_usol, atol=tol_usol,dtmax=dtmax,
        )
    end
    opt1 = find_zero(fun, (lb, ub), Roots.Brent(); atol=tol_A0, rtol=tol_A0)
    @show optmin1 = opt1
    sol = du1solve(
        O1du1,
        optmin1;
        rhomin=rhomin,
        rhomax=rhomax,
        rtol=tol_usol,
        atol=tol_usol,
        d=d,
        n=n,
        η=η,
        dtmax=dtmax,
        method=RadauIIA5(),
    )
    return USolution(sol.sol; A0=optmin1, rhomin=rhomin, rhomax=rhomax, d=d, n=n, η=η)
end

function findeigenval(
    Eigenfun; Usol, lambdalb::T, lambdaub::T, tol_eigensol, tol_λ, method=Brent()
) where {T}
    u1fun = Interpolations.interpolate((Usol.t,), Usol.u2, Gridded(Linear()))
    u2fun = Interpolations.interpolate((Usol.t,), Usol.u3, Gridded(Linear()))
    function fun(para)
        @show para
        loss = LossEigensolve(
            Eigenfun;
            u1fun=u1fun,
            u2fun=u2fun,
            rhomin=Usol.rhomin,
            rhomax=Usol.rhomax,
            A0=T(Usol.A0),
            d=Usol.d,
            n=Usol.n,
            η=Usol.η,
            λ=para,
            atol=tol_eigensol,
            rtol=tol_eigensol,
            method=RadauIIA5(),
        )
        return @show loss
    end
    lb = lambdalb
    ub = lambdaub
    opt1 = Optim.optimize(fun, lb, ub, method; rel_tol=T(tol_λ), abs_tol=T(tol_λ))
    @show optmin1 = opt1.minimizer
    return USolution(
        Usol.sol;
        A0=Usol.A0,
        rhomin=Usol.rhomin,
        rhomax=Usol.rhomax,
        d=Usol.d,
        n=Usol.n,
        η=Usol.η,
        λ=optmin1,
    )
end

function findeigenval2(
    Eigenfun; Usol, lambdalb::T, lambdaub::T, tol_eigensol, tol_λ, method=Brent()
) where {T}
    u1fun = Interpolations.interpolate((Usol.t,), Usol.u2, Gridded(Linear()))
    u2fun = Interpolations.interpolate((Usol.t,), Usol.u3, Gridded(Linear()))
    function fun(para)
        @show para
        loss = LossEigensolve2(
            Eigenfun;
            u1fun=u1fun,
            u2fun=u2fun,
            rhomin=Usol.rhomin,
            rhomax=Usol.rhomax,
            A0=-Usol.A0,
            d=Usol.d,
            n=Usol.n,
            η=Usol.η,
            λ=para,
            atol=tol_eigensol,
            rtol=tol_eigensol,
            method=RadauIIA5(),
        )
        return @show loss
    end
    lb = lambdalb
    ub = lambdaub
    # opt1 = Optim.optimize(fun, lb, ub, method; rel_tol=T(1e-15), abs_tol=T(1e-15))
    # @show optmin1 = opt1.minimizer
    optmin1 = find_zero(fun, (lb, ub), Bisection(); atol=tol_λ, rtol=tol_λ)
    return USolution(
        Usol.sol;
        A0=Usol.A0,
        rhomin=Usol.rhomin,
        rhomax=Usol.rhomax,
        d=Usol.d,
        n=Usol.n,
        η=Usol.η,
        λ=optmin1,
    )
end

function findeigenval3(
    Eigenfun; Usol, lambdalb::T, lambdaub::T, tol_eigensol, tol_λ, method=Brent()
) where {T}
    u1fun = Interpolations.interpolate((Usol.t,), Usol.u2, Gridded(Linear()))
    u2fun = Interpolations.interpolate((Usol.t,), Usol.u3, Gridded(Linear()))
    function fun(para)
        @show para
        loss = LossEigensolve2(
            Eigenfun;
            u1fun=u1fun,
            u2fun=u2fun,
            rhomin=Usol.rhomin*100,
            rhomax=Usol.rhomax,
            A0=-Usol.A0*10^10,
            d=Usol.d,
            n=Usol.n,
            η=Usol.η,
            λ=para,
            atol=tol_eigensol,
            rtol=tol_eigensol,
            method=RadauIIA5(),
        )
        return @show loss
    end
    myT=eltype([Usol.A0,Usol.d,Usol.n,Usol.η])
    lb = lambdalb|>myT
    ub = lambdaub|>myT
    # opt1 = Optim.optimize(fun, lb, ub, method; rel_tol=T(1e-15), abs_tol=T(1e-15))
    # @show optmin1 = opt1.minimizer
    optmin1 = find_zero(fun, (lb, ub), Bisection(); atol=tol_λ, rtol=tol_λ)
    tspan = (Usol.rhomax, Usol.rhomin*100)
    p = (Usol.d, Usol.n, Usol.η, optmin1)
    u0 = iniEigen(Usol.A0*10^10, Usol.rhomax, Usol.d, Usol.η, optmin1)
    f = (du, u, p, rho) -> Eigenfun(du, u, p, rho, u1fun, u2fun)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(
        prob, method; abstol=tol_eigensol, reltol=tol_eigensol, maxiters=10^7
    )
    return USolution(
        sol;
        A0=Usol.A0,
        rhomin=Usol.rhomin*100,
        rhomax=Usol.rhomax,
        d=Usol.d,
        n=Usol.n,
        η=Usol.η,
        λ=optmin1,
    )
end




function findeigenval2_y2(
    Eigenfun; Usol, lambdalb::T, lambdaub::T, tol_eigensol, tol_λ, method=Brent()
) where {T}
    myT=typeof(Usol.A0)
    u1fun = Interpolations.interpolate((Usol.t,), Usol.u2, Gridded(Linear()))
    u2fun = Interpolations.interpolate((Usol.t,), Usol.u3, Gridded(Linear()))
    function fun(para)
        @show para
        loss = LossEigensolve2(
            Eigenfun;
            u1fun=u1fun,
            u2fun=u2fun,
            rhomin=Usol.rhomin*10,
            rhomax=Usol.rhomax,
            A0=-Usol.A0,
            d=Usol.d,
            n=Usol.n,
            η=Usol.η,
            λ=para,
            atol=tol_eigensol,
            rtol=tol_eigensol,
            method=RadauIIA5(),
        )
        return @show loss
    end
    lb = lambdalb|>myT
    ub = lambdaub|>myT
    # opt1 = Optim.optimize(fun, lb, ub, method; rel_tol=T(1e-15), abs_tol=T(1e-15))
    # @show optmin1 = opt1.minimizer
    optmin1 = find_zero(fun, (lb, ub), Bisection(); atol=tol_λ, rtol=tol_λ)
    return optmin1
end




function Udatasaver!(udata::Dict; d, n, lb_A0, ub_A0, lb_λ, ub_λ)
    sol = findUmin(;
        lb=myT(lb_A0),
        ub=myT(ub_A0),
        rhomin=myT(10^-8),
        rhomax=myT(10),
        d=myT(d),
        n=myT(n),
        η=myT(0),
        atol=1e-20,
        rtol=1e-20,
    )
    usol = findeigenval(
        Eigenfun;
        Usol=sol,
        lambdalb=myT(lb_λ),
        lambdaub=myT(ub_λ),
        rtol=1e-16,
        atol=1e-16,
        method=Brent(),
    )
    return udata[d, n] = usol
end


function Udatasaver2!(
    udata::Dict;
    d,
    n,
    lb_A0,
    ub_A0,
    lb_λ,
    ub_λ,
    rhomin_pre=1//100,
    tol_usol,
    tol_A0,
    tol_λ,
    tol_eigensol,
)
    myT = eltype([d, n, lb_A0, ub_A0, lb_λ, ub_λ])
    sol = findUmin(;
        lb=myT(lb_A0),
        ub=myT(ub_A0),
        rhomin=myT(rhomin_pre),
        rhomax=myT(10),
        d=myT(d),
        n=myT(n),
        η=myT(0),
        tol_usol=myT(tol_usol),
        tol_A0=myT(tol_A0),
    )
    usol = findeigenval(
        Eigenfun;
        Usol=sol,
        lambdalb=myT(lb_λ),
        lambdaub=myT(ub_λ),
        method=Brent(),
        tol_eigensol=myT(tol_eigensol),
        tol_λ=myT(tol_λ),
    )
    return udata[d, n] = usol
end



function Udatasaver3!(
    udata::Dict;
    d,
    n,
    lb_A0,
    ub_A0,
    lb_λ,
    ub_λ,
    η=0,
    rhomin_pre=1//1000000,
    tol_usol,
    tol_A0,
    tol_λ,
    tol_eigensol,
    dtmax=0.005,
)
    myT = eltype([d, n, lb_A0, ub_A0, lb_λ, ub_λ])
    # sol_pre = findUmin(;
    #     lb=myT(lb_A0),
    #     ub=myT(ub_A0),
    #     rhomin=myT(rhomin_pre),
    #     rhomax=myT(10),
    #     d=myT(d),
    #     n=myT(n),
    #     η=myT(0),
    #     tol_usol=myT(tol_usol),
    #     tol_A0=myT(tol_A0*10^10),
    # )
    sol = findUmin2(;
        lb=myT(lb_A0),
        ub=myT(ub_A0),
        rhomin=myT(rhomin_pre),
        rhomax=myT(100),
        d=myT(d),
        n=myT(n),
        η=myT(η),
        tol_usol=myT(tol_usol),
        tol_A0=myT(tol_A0),
        dtmax=dtmax,
    )
    usol = findeigenval(
        Eigenfun;
        Usol=sol,
        lambdalb=myT(lb_λ),
        lambdaub=myT(ub_λ),
        method=Brent(),
        tol_eigensol=myT(tol_eigensol),
        tol_λ=myT(tol_λ),
    )
    d1p=rationalize(Float32(d))
    n1p=rationalize(Float32(n))
    return udata[d1p, n1p] = usol
end

function minisave!(;udic_mini::Dict,udic_large::Dict)
    for (d, n) in keys(udic_large)
        udic_mini[d,n]=USolution_mini(usol_WF[d, n])
    end
    return nothing
end

function etasave(path,usol)
    rho0=find_zero(x->usol.sol(x)[2], (0.00001, 1.0))
    upp=usol.sol(rho0)[3]
    uppp=usol.sol(rho0,Val{1})[3]
    v=[usol.d,usol.n,usol.η,rho0,upp,uppp]
    writedlm(path,v)
end
