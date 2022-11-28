include("load.jl")
include("struct.jl")
include("ini.jl")
include("regulator.jl")
include("flow.jl")
include("solver.jl")
myT = Double64
myT2 = Float128
setprecision(BigFloat, 128)
@time O1sol = O1du1solve(
    O1du1,
    myT(val_A0);
    rhomin=myT(0.00001),
    rho0=myT(10),
    rtol=1e-20,
    atol=1e-20,
    method=RadauIIA5(),
)

@time O1du0sol = du0solve(
    O1du0, myT(val_A0); rtol=1e-16, atol=1e-16, method=RadauIIA5(), d=myT(3)
)

@time O1du0sol = du0solve(
    O1du0, myT(val_A0); rtol=1e-16, atol=1e-16, method=RadauIIA5(), d=myT(3)
)

val_A0 = find_zero(getA0, (myT(84.1 / 3), myT(84.2 / 3)), Bisection(); rtol=1e-4, atol=1e-4)
val_A02 = find_zero(
    getA02, (myT(val_A0 - 1e-9), myT(val_A0 + 1e-9)), Bisection(); rtol=1e-2, atol=1e-2
)

1 / val_lambda
interpolate(O1sol.t, O1sol.u3, BSplineOrder(1))

itp1 = Interpolations.interpolate((O1du0sol.t,), O1du0sol.u2, Gridded(Linear()))
itp2 = Interpolations.interpolate((O1du0sol.t,), O1du0sol.u3, Gridded(Linear()))

eigensol = Eigensolve(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    rhomin=myT(0.00001),
    rho0=myT(10),
    A0=myT(-170),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    λ=myT(1.6),
    rtol=1e-14,
    atol=1e-14,
    method=RadauIIA5(),
)

opteigensol = getlambda0(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    lambdalb=myT(-6.0),
    deltalambda=myT(0.5),
    rhomin=myT(0.00001),
    rho0=myT(10),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    rtol=1e-12,
    atol=1e-12,
    method=RadauIIA5(),
)
1
minimum(opteigensol)
1 / Optim.minimizer(opteigensol)[2]
Optim.minimizer(opteigensol)

plot(eigensol.t[1:100], eigensol.u1[1:100])
plot(eigensol.t[1:200], eigensol.u2[1:200])

plot!(O1sol.t[1:30], O1sol.u2[1:30])
plot(O1du0sol.t[1:300], O1du0sol.u3[1:300])

plot(eigensol.t[1:100], eigensol.u2[1:100])

plot(O1sol.t[500:3000], geterror(O1sol)[500:3000])

plot(O1sol.t[1:40], O1sol.u2[1:40])

plot!(O1sol2.t[1:10], O1sol2.u3[1:10])

writedlm("data/du1_d=3N=1.dat", [O1sol.t O1sol.u1 O1sol.u2 O1sol.u3])

writedlm("data/du0_d=3N=1.dat", [O1du0sol.t O1du0sol.u1 O1du0sol.u2 O1du0sol.u3])

1

O1sol.t

plot!(taylor, 0, 0.005)

function taylor(ρ)
    return -0.18616963580414798662258800563977391 +
           4.8678407813083455712411687656558631 * ρ +
           33.31749708986721867035895808518545 * ρ^2 +
           182.3300122955358880765231767203499 * ρ^3 +
           572.191967694137673274069617988882 * ρ^4 -
           1476.7212074253873855829730305918 * ρ^5 -
           30574.074925225317423274882460331 * ρ^6 -
           119358.99387510873714047993885044 * ρ^7 +
           1102114.758259733917713234579554 * ρ^8 +
           18840284.21530451336727325919163 * ρ^9 +
           111951150.6595689862048869282864 * ρ^10 +
           82947590.46448809032702176621 * ρ^11
end
