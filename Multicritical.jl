#W-F fix point
val_A0_mul1 = find_zero(
    x -> getA02(x; d=myT(2.4), rhomax=myT(100)),
    (myT(2000), myT(4000)),
    Bisection();
    rtol=1000,
    atol=1000,
)
val_A0_mul1_1 = find_zero(
    x -> getA02(x; d=myT2(2.4), rhomax=myT2(100)),
    (myT2(val_A0_mul1 - 1e-10), myT2(val_A0_mul1 + 1e-10)),
    Bisection();
    rtol=300,
    atol=300,
)

sol_O1_d24_mul1 = du1solve(
    O1du1,
    myT2(val_A0_mul1_1);
    rhomin=myT2(1e-6),
    rhomax=myT2(100),
    rtol=1e-22,
    atol=1e-22,
    d=myT2(2.4),
    dtmax=0.005,
    method=RadauIIA5(),
)
plot(sol_O1_d24_mul1.t[1:8000], sol_O1_d24_mul1.u1[1:8000])

itp1_O1_d24_mul1 = Interpolations.interpolate(
    (sol_O1_d24_mul1.t,), sol_O1_d24_mul1.u2, Gridded(Linear())
)
itp2_O1_d24_mul1 = Interpolations.interpolate(
    (sol_O1_d24_mul1.t,), sol_O1_d24_mul1.u3, Gridded(Linear())
)
#
eigensol_lambda1_O1_d24_mul1 = Eigensolve(
    Eigenfun;
    u1fun=itp1_O1_d24_mul1,
    u2fun=itp2_O1_d24_mul1,
    rhomin=myT2(1e-6),
    rhomax=myT2(100),
    A0=myT2(val_A0_mul1_1),
    d=myT2(2.4),
    n=myT2(1),
    η=myT2(0),
    λ=myT2(opt_lambda1_O1_d24_mul1),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)
plot(eigensol_O1_d24_mul1.t[1:40], eigensol_O1_d24_mul1.u1[1:40])

opt_lambda1_O1_d24_mul1 = getlambda0(
    Eigenfun;
    u1fun=itp1_O1_d24_mul1,
    u2fun=itp2_O1_d24_mul1,
    lambdalb=myT2(0.9),
    lambdaub=myT2(1.2),
    rhomin=myT2(1e-6),
    rhomax=myT2(100),
    A0=myT2(val_A0_mul1_1),
    d=myT2(2.4),
    n=myT2(1),
    η=myT2(0),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)
#
eigensol_lambda2_O1_d24_mul1 = Eigensolve(
    Eigenfun;
    u1fun=itp1_O1_d24_mul1,
    u2fun=itp2_O1_d24_mul1,
    rhomin=myT(1e-6),
    rhomax=myT(100),
    A0=myT(val_A0_mul1_1),
    d=myT(2.4),
    n=myT(1),
    η=myT(0),
    λ=myT(opt_lambda2_O1_d24_mul1),
    rtol=1e-12,
    atol=1e-12,
    method=RadauIIA5(),
)
plot(eigensol_lambda2_O1_d24_mul1.t[1:40], eigensol_lambda2_O1_d24_mul1.u2[1:40])

opt_lambda2_O1_d24_mul1 = getlambda0(
    Eigenfun;
    u1fun=itp1_O1_d24_mul1,
    u2fun=itp2_O1_d24_mul1,
    lambdalb=myT(-1.0),
    lambdaub=myT(-0.9),
    rhomin=myT(1e-6),
    rhomax=myT(100),
    A0=myT(val_A0_mul1_1),
    d=myT(2.4),
    n=myT(1),
    η=myT(0),
    rtol=1e-14,
    atol=1e-14,
    method=RadauIIA5(),
)

#
eigensol_lambda3_O1_d24_mul1 = Eigensolve(
    Eigenfun;
    u1fun=itp1_O1_d24_mul1,
    u2fun=itp2_O1_d24_mul1,
    rhomin=myT(1e-6),
    rhomax=myT(100),
    A0=myT(val_A0_mul1_1),
    d=myT(2.4),
    n=myT(1),
    η=myT(0),
    λ=myT(opt_lambda3_O1_d24_mul1),
    rtol=1e-12,
    atol=1e-12,
    method=RadauIIA5(),
)
plot(
    eigensol_lambda3_O1_d24_mul1.t[1:600],
    eigensol_lambda3_O1_d24_mul1.u1[1:600] ./ eigensol_lambda3_O1_d24_mul1.u1[1],
)

opt_lambda3_O1_d24_mul1 = getlambda0(
    Eigenfun;
    u1fun=itp1_O1_d24_mul1,
    u2fun=itp2_O1_d24_mul1,
    lambdalb=myT(-3.5),
    lambdaub=myT(-3.0),
    rhomin=myT(1e-6),
    rhomax=myT(100),
    A0=myT(val_A0_mul1_1),
    d=myT(2.4),
    n=myT(1),
    η=myT(0),
    rtol=1e-14,
    atol=1e-14,
    method=RadauIIA5(),
)

#
eigensol_lambda4_O1_d24_mul1 = Eigensolve(
    Eigenfun;
    u1fun=itp1_O1_d24_mul1,
    u2fun=itp2_O1_d24_mul1,
    rhomin=myT(1e-6),
    rhomax=myT(100),
    A0=myT(val_A0_mul1_1),
    d=myT(2.4),
    n=myT(1),
    η=myT(0),
    λ=myT(-4.0),
    rtol=1e-12,
    atol=1e-12,
    method=RadauIIA5(),
)
plot(eigensol_lambda4_O1_d24_mul1.t[1:40], eigensol_lambda4_O1_d24_mul1.u2[1:40])

opt_lambda4_O1_d24_mul1 = getlambda0(
    Eigenfun;
    u1fun=itp1_O1_d24_mul1,
    u2fun=itp2_O1_d24_mul1,
    lambdalb=myT(-3.5),
    lambdaub=myT(-3.0),
    rhomin=myT(1e-6),
    rhomax=myT(100),
    A0=myT(val_A0_mul1_1),
    d=myT(2.4),
    n=myT(1),
    η=myT(0),
    rtol=1e-14,
    atol=1e-14,
    method=RadauIIA5(),
)

#3-critical fix point
val_A0_mul2 = find_zero(
    x -> getA02(x; d=myT(2.4), rhomax=myT(100)),
    (myT(0.1), myT(2.0)),
    Bisection();
    rtol=1000,
    atol=1000,
)
val_A0_mul2_1 = find_zero(
    x -> getA02(x; d=myT2(2.4), rhomax=myT2(100)),
    (myT2(val_A0_mul2 - 1e-10), myT2(val_A0_mul2 + 1e-10)),
    Bisection();
    rtol=10000000,
    atol=10000000,
)

sol_O1_d24_mul2 = du1solve(
    O1du1,
    myT2(val_A0_mul2_1);
    rhomin=myT2(1e-6),
    rhomax=myT2(100),
    rtol=1e-22,
    atol=1e-22,
    d=myT2(2.4),
    dtmax=0.005,
    method=RadauIIA5(),
)
plot(sol_O1_d24_mul2.t[1:13000], sol_O1_d24_mul2.u1[1:13000])

itp1_O1_d24_mul2 = Interpolations.interpolate(
    (sol_O1_d24_mul2.t,), sol_O1_d24_mul2.u2, Gridded(Linear())
)
itp2_O1_d24_mul2 = Interpolations.interpolate(
    (sol_O1_d24_mul2.t,), sol_O1_d24_mul2.u3, Gridded(Linear())
)

#
eigensol_lambda1_O1_d24_mul2 = Eigensolve(
    Eigenfun;
    u1fun=itp1_O1_d24_mul2,
    u2fun=itp2_O1_d24_mul2,
    rhomin=myT2(1e-6),
    rhomax=myT2(100),
    A0=myT2(val_A0_mul2_1),
    d=myT2(2.4),
    n=myT2(1),
    η=myT2(0),
    λ=myT2(opt_lambda1_O1_d24_mul2),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)
plot(eigensol_lambda1_O1_d24_mul2.t[1:40], eigensol_lambda1_O1_d24_mul2.u1[1:40])

opt_lambda1_O1_d24_mul2 = getlambda0(
    Eigenfun;
    u1fun=itp1_O1_d24_mul2,
    u2fun=itp2_O1_d24_mul2,
    lambdalb=myT2(0.9),
    lambdaub=myT2(1.2),
    rhomin=myT2(1e-6),
    rhomax=myT2(100),
    A0=myT2(val_A0_mul2_1),
    d=myT2(2.4),
    n=myT2(1),
    η=myT2(0),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)

#4-critical fix point
val_A0_mul3 = find_zero(
    x -> getA02(x; d=myT(2.4), rhomax=myT(100)),
    (myT(0.005), myT(0.1)),
    Bisection();
    rtol=1000,
    atol=1000,
)
val_A0_mul3_1 = find_zero(
    x -> getA02(x; d=myT2(2.4), rhomax=myT2(100)),
    (myT2(val_A0_mul3 - 1e-10), myT2(val_A0_mul3 + 1e-10)),
    Bisection();
    rtol=10000000,
    atol=10000000,
)

sol_O1_d24_mul3 = du1solve(
    O1du1,
    myT2(val_A0_mul3_1);
    rhomin=myT2(1e-6),
    rhomax=myT2(100),
    rtol=1e-22,
    atol=1e-22,
    d=myT2(2.4),
    dtmax=0.005,
    method=RadauIIA5(),
)
plot(sol_O1_d24_mul3.t[1:8000], sol_O1_d24_mul3.u1[1:8000])

itp1_O1_d24_mul3 = Interpolations.interpolate(
    (sol_O1_d24_mul3.t,), sol_O1_d24_mul3.u2, Gridded(Linear())
)
itp2_O1_d24_mul3 = Interpolations.interpolate(
    (sol_O1_d24_mul3.t,), sol_O1_d24_mul3.u3, Gridded(Linear())
)

#5-critical fix point
val_A0_mul4 = find_zero(
    x -> getA02(x; d=myT(2.4), rhomax=myT(100)),
    (myT(0.0001), myT(0.005)),
    Bisection();
    rtol=1000,
    atol=1000,
)
val_A0_mul4_1 = find_zero(
    x -> getA02(x; d=myT2(2.4), rhomax=myT2(100)),
    (myT2(val_A0_mul4 - 1e-15), myT2(val_A0_mul4 + 1e-15)),
    Bisection();
    rtol=2 * 10^6,
    atol=2 * 10^6,
)

sol_O1_d24_mul4 = du1solve(
    O1du1,
    myT2(val_A0_mul4_1);
    rhomin=myT2(1e-6),
    rhomax=myT2(100),
    rtol=1e-22,
    atol=1e-22,
    d=myT2(2.4),
    dtmax=0.005,
    method=RadauIIA5(),
)
plot(sol_O1_d24_mul4.t[1:4000], sol_O1_d24_mul4.u1[1:4000])

itp1_O1_d24_mul4 = Interpolations.interpolate(
    (sol_O1_d24_mul4.t,), sol_O1_d24_mul4.u2, Gridded(Linear())
)
itp2_O1_d24_mul4 = Interpolations.interpolate(
    (sol_O1_d24_mul4.t,), sol_O1_d24_mul4.u3, Gridded(Linear())
)

writedlm(
    "data/O1/d=2.4/sol_O1_d24_mul1.dat",
    [sol_O1_d24_mul1.t sol_O1_d24_mul1.u1 sol_O1_d24_mul1.u2 sol_O1_d24_mul1.u3],
)
writedlm(
    "data/O1/d=2.4/sol_O1_d24_mul2.dat",
    [sol_O1_d24_mul2.t sol_O1_d24_mul2.u1 sol_O1_d24_mul2.u2 sol_O1_d24_mul2.u3],
)
writedlm(
    "data/O1/d=2.4/sol_O1_d24_mul3.dat",
    [sol_O1_d24_mul3.t sol_O1_d24_mul3.u1 sol_O1_d24_mul3.u2 sol_O1_d24_mul3.u3],
)
writedlm(
    "data/O1/d=2.4/sol_O1_d24_mul4.dat",
    [sol_O1_d24_mul4.t sol_O1_d24_mul4.u1 sol_O1_d24_mul4.u2 sol_O1_d24_mul4.u3],
)

writedlm("data/O1/d=2.4/val_A0_mul1.dat", val_A0_mul1_1)
writedlm("data/O1/d=2.4/val_A0_mul2.dat", val_A0_mul2_1)
writedlm("data/O1/d=2.4/val_A0_mul3.dat", val_A0_mul3_1)
writedlm("data/O1/d=2.4/val_A0_mul4.dat", val_A0_mul4_1)


plot(eigensol_lambda1_O1_d24_mul1.t[1:100], eigensol_lambda1_O1_d24_mul1.u1[1:100]./eigensol_lambda1_O1_d24_mul1.u1[1])
plot!(eigensol_lambda2_O1_d24_mul1.t[1:380], eigensol_lambda2_O1_d24_mul1.u1[1:380]./eigensol_lambda2_O1_d24_mul1.u1[1])
plot!(
    eigensol_lambda3_O1_d24_mul1.t[1:650],
    eigensol_lambda3_O1_d24_mul1.u1[1:650] ./ eigensol_lambda3_O1_d24_mul1.u1[1],
)
plot(
    eigensol_lambda4_O1_d24_mul1.t[1:380],
    eigensol_lambda4_O1_d24_mul1.u1[1:380] ./ eigensol_lambda4_O1_d24_mul1.u1[1],
)
