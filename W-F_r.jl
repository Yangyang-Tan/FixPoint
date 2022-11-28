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
usol_r4 = Dict()
usol_r_MWH = Dict()
# n=1
usol_r1 = Dict()
myT = Float64
Udatasaver3!_r(
    usol_r4;
    d=2.0,
    n=1.0,
    η=0.3,
    lb_A0=myT(0.5),
    ub_A0=myT(1000.0),
    lb_λ=0.5,
    ub_λ=1.8,
    rhomin_pre=1//100,
    tol_usol=1e-15,
    tol_A0=1e-20,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    nr=2,
)

plot(x -> usol_r4[(2//1, 1//1)].sol(x)[2], 0.01, 0.2; label="n=1")

etafun_r(usol_r4[(2//1, 1//1)], x->exp(-x)/x, x->-(((1+x)*exp(-x))/x^2), x->((2+2*x+x^2)*exp(-x))/x^3)


etafun_r(usol_r4[(2//1, 1//1)], r_exp(2), rp_exp(2), rpp_exp(2))


plot(x -> usol_r1[(3//1, 1//1)].sol(x)[2], 0.0001, 0.1; label="n=1")
plot!(x -> usol_r2[(3//1, 1//1)].sol(x)[2], 0.0001, 0.1; label="n=2")
plot!(x -> usol_r4[(3//1, 1//1)].sol(x)[2], 0.0001, 0.1; label="n=4")
plot!(x -> usol_WF[(3.0, 1.0)].sol(x)[2], 0.0001, 0.1; label="Opt")

plot(
    1 ./ [
        usol_r1[(3//1, 1//1)].λ,
        usol_r2[(3//1, 1//1)].λ,
        usol_r4[(3//1, 1//1)].λ,
        usol_WF[(3.0, 1.0)].λ,
    ];
    seriestype=:bar,
    label="λ",
    xticks=(1:4, ["n=1", "n=2", "n=4", "Opt"]),
)

testf1 = lpfunSp(myT(3.0), 2.0, 0.0)
testf2 = lpfun(3.0, 2.0, 0.0)

xrang = (0.1:0.1:100)
@time testf2.(xrang)
@time testf1.(xrang)

plot(x -> testf1(x), 0.0, 10.0; label="n=2")
plot!(x -> testf2(x), 0.1, 10.0; label="n=2")
