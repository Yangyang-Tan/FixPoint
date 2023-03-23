tempsol = du1solve(
    effective_r,
    myT(10.8);
    rhomin=myT(0.001),
    rhomax=myT(100.0),
    d=myT(3.0),
    n=myT(1),
    η=myT(0.0),
    atol=1e-15,
    rtol=1e-15,
    dtmax=0.005,
    method=RadauIIA5(),
    inifun=ini_CS,
)
usol_r2 = Dict()
usol_r_MWH = Dict()
# n=1
usol_r1 = Dict()
myT = Double64

Udatasaver3!_r(
    usol_r2;
    d=2.0,
    n=1.7,
    η=0.15,
    lb_A0=myT(100.0),
    ub_A0=myT(1999.75),
    lb_λ=0.1,
    ub_λ=1.8,
    rhomin_pre=1//100,
    tol_usol=1e-12,
    tol_A0=1e-12,
    tol_eigensol=1e-12,
    tol_λ=1e-12,
    nr=2,
)

plot(x -> usol_r2[(2//1, 17//10)].sol(x)[2], 0.1, 0.6; label="n=1")

myT = Float64

for j in 1:3
    Udatasaver3!_r(
        usol_r2;
        d=2.0,
        n=1.7,
        η=etafun_r(usol_r2[(2//1, 17//10)], r_exp(2), rp_exp(2), rpp_exp(2)),
        lb_A0=myT(100.0),
        ub_A0=myT(2200.75),
        lb_λ=0.1,
        ub_λ=1.0,
        rhomin_pre=1//100,
        tol_usol=1e-12,
        tol_A0=1e-5,
        tol_eigensol=1e-12,
        tol_λ=1e-10,
        nr=2,
    )
end

tempsol = du1solve(
    effective_rfun(myT(2.0), myT(2.0), myT(0.04588436529772546)),
    myT(1.7039554750420893e-11);
    rhomin=myT(0.1),
    rhomax=myT(20.0),
    d=myT(2.0),
    n=myT(1.9),
    η=myT(0.04588436529772546),
    atol=1e-21,
    rtol=1e-21,
    dtmax=0.005,
    method=RadauIIA5(),
    inifun=ini_CS,
)
tempsol.d = 2.0
plot(x -> tempsol.sol(x)[2], 0.1, 1.7; label="n=1")
tempsol.sol(0.2)[2]
etafun_r(usol_r2[(2//1, 10//10)], r_exp(2), rp_exp(2), rpp_exp(2))
etafun_r(tempsol, r_exp(2), rp_exp(2), rpp_exp(2))

usol_r2[(2//1, 19//10)].λ
plot(
    2.0 .- [1.0, 1.5, 1.6, 1.7, 1.8, 1.9],
    [
        usol_r2[(2//1, 10//10)].λ,
        usol_r2[(2//1, 15//10)].λ,
        usol_r2[(2//1, 16//10)].λ,
        usol_r2[(2//1, 17//10)].λ,
        usol_r2[(2//1, 18//10)].λ,
        usol_r2[(2//1, 19//10)].λ,
    ];
    label="",
    seriestype=:scatter,
    xaxis=:log,
    yaxis=:log,
    xlabel="2-N",
    ylabel="1/ν"
)

plot(
    2.0 .- [1.0, 1.5, 1.6, 1.7, 1.8, 1.9],
    [
        usol_r2[(2//1, 10//10)].η,
        usol_r2[(2//1, 15//10)].η,
        usol_r2[(2//1, 16//10)].η,
        usol_r2[(2//1, 17//10)].η,
        usol_r2[(2//1, 18//10)].η,
        usol_r2[(2//1, 19//10)].η,
    ];
    label="",
    seriestype=:scatter,
    xaxis=:log,
    yaxis=:log,
    xlabel="2-N",
    ylabel="η",
)
