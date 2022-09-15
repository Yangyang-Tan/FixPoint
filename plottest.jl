plot(O1du0sol.t[1:1000], O1du0sol.u2[1:1000])
plot(O1sol.t[1:600], O1sol.u3[1:600])
plot(O1sol.t[1:13000], O1sol.u1[1:13000])

Nderivative(O1sol.t[1:5000], O1sol.u3[1:5000], 2*1e-5)

@time O1sol = du1solve(
    O1du1,
    myT2(val_A02+1e-8);
    rhomin=myT2(0.000000001),
    rho0=myT2(10),
    rtol=1e-22,
    atol=1e-22,
    d=myT2(2.4),
    dtmax=0.005,
    method=RadauIIA5(),
)
val_A02 = find_zero(x->getA02(x,d=myT2(2.4)), (myT2(1.0), myT2(1.1)), Bisection(), rtol=1e-14, atol=1e-14)

val_A021 = find_zero(x->getA02(x,d=myT2(2.4)), (myT2(val_A02-1e-12), myT2(val_A02+1e-12)), Bisection(), rtol=100, atol=100)

O1sol.t
