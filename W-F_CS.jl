usol_CS = Dict()
usol_CS_MWH = Dict()
# n=1
Udatasaver3!_CS(
    usol_CS;
    d=2.0,
    n=0.8,
    η=0.07,
    lb_A0=myT(0.53890),
    ub_A0=myT(0.538905),
    lb_λ=0.1,
    ub_λ=1.0,
    rhomin_pre=1//10,
    tol_usol=1e-18,
    tol_A0=1e-10,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
)
8.572500445465838e-7
for i in 1:1
    n = rationalize(0.9)
    d = rationalize(2.0)
    Udatasaver3!_CS(
        usol_CS;
        d=d,
        n=n,
        η=etafun_CS(usol_CS[d, n], d),
        lb_A0=0.01 * usol_CS[d, n].A0,
        ub_A0=10.2 * usol_CS[d, n].A0,
        lb_λ=0.05,
        ub_λ=1.5,
        rhomin_pre=1//100,
        tol_usol=1e-18,
        tol_A0=1e-18,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end

temsoleigen = Eigensolve(
    Eigenfun_CS;
    u1fun=tempu1fun,
    u2fun=tempu2fun,
    rhomin=myT(0.001),
    rhomax=myT(50),
    A0=myT(42.06),
    d=myT(3.0),
    n=myT(1.0),
    η=myT(0.2),
    λ=myT(-3.0),
    atol=1e-16,
    rtol=1e-16,
    method=RadauIIA5(),
)

plot(x -> temsoleigen.sol(x)[2], 0.0001, 0.1)

@time tempsol = du1solve(
    effective_CS,
    myT(10.8);
    rhomin=myT(0.01),
    rhomax=myT(100.0),
    d=myT(2.0),
    n=myT(1),
    η=myT(0.1),
    atol=1e-21,
    rtol=1e-21,
    dtmax=0.005,
    method=RadauIIA5(),
    inifun=ini_CS,
)

plot(x -> tempsol.sol(x, Val{1})[3], 0.01, 1)
plot(x -> tempu2fun(x), 0.0001, 0.1)
plot(tempsol.t, tempsol.u3)
ini_CS(1.0, 100.0)
inidu0(1.0, 100.0)
tempdu2 = randn(2)
tempdu1 = randn(3)
effective_CS2(
    tempdu2,
    ini_CS2(2.0, 99.0),
    [3.0, 1.0, 0.0, (2) * gamma(1 - 3 / 2) / (2 * (4 * π)^(3 / 2))],
    99.0,
)
O1du1(
    tempdu1,
    inidu0(2.0, 100.0),
    [3.0, 1.0, 0.0, (2 + 3) / (pi^(3 / 2) * (3 * (2 + 3) * gamma(3 / 2) * 2^(3 - 1)))],
    100.0,
)

ini_CS2(2.0, 100.0)
inidu0(2.0, 100.0)
tempdu2
tempdu1
[3.0, 1.0, 0.0, (2) * gamma(1 - 3 / 2) / (2 * (4 * π)^(3 / 2))]
ini_CS2(2.0, 99.0)
(2) * gamma(1 - 3 / 2) / (2 * (4 * π)^(3 / 2))
(2 + 3) / (pi^(3 / 2) * (3 * (2 + 3) * gamma(3 / 2) * 2^(3 - 1)))
