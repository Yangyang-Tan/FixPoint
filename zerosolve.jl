tempsol=let d=myT(30//20), n=myT(1//1), η=myT(0.0), A0=0.0, rhomin=myT(0.01), dtmax=0.005
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    # c1 = (2-η)*gamma(2-d/2)/(2*(4*π)^(d/2))
    p = [d, n, η, c1]
    f = ODEFunction(O1du1)
    tspan = (rhomin,10*usol_MWH_d2n2_eta[rationalize(0.1)].rhomax)
    # differential_vars = [true, true, false]
    # u0 = iniO1du1(A0, rhomax)
    # u0 = inifun(A0, rhomax; d=d, n=n, η=η)
    σ=myT(-0.35)
    u0 = [(c1/d)/(1+σ),σ,σ]
    prob = ODEProblem(f, u0, tspan, p)
    solve(
        prob, RadauIIA5(); abstol=1e-6, reltol=1e-6, maxiters=10^7, dtmin=1e-27
    )
    # return USolution(sol)
end
tempsol.sol

plot(x->tempsol(x)[2],0.01,0.1)
usol_MWH_d2n2_eta[rationalize(0.058)].t[1]
usol_MWH_d2n2_eta[rationalize(0.064)].u1[1]
usol_MWH_d2n2_eta[rationalize(0.064)].u2[1]
