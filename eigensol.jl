opteigensol1 = getlambda0(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    lambdalb=myT(1.5),
    deltalambda=myT(0.1),
    rhomin=myT(0.00001),
    rho0=myT(10),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    rtol=1e-12,
    atol=1e-12,
    method=RadauIIA5(),
)

itp1d24 = Interpolations.interpolate((O1sol.t,), O1sol.u2, Gridded(Linear()))
itp2d24 = Interpolations.interpolate((O1sol.t,), O1sol.u3, Gridded(Linear()))

opteigensol1d24 = getlambda0(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    lambdalb=myT2(1.5),
    deltalambda=myT2(0.1),
    rhomin=myT2(0.00001),
    rho0=myT2(10),
    d=myT2(2.4),
    n=myT2(1),
    η=myT2(0),
    rtol=1e-12,
    atol=1e-12,
    method=RadauIIA5(),
)

eigensol1d24 = Eigensolve(
    Eigenfun;
    u1fun=itp1d24,
    u2fun=itp2d24,
    rhomin=myT2(1e-8),
    rho0=myT2(10),
    A0=myT2(10.0),
    d=myT2(3.0),
    n=myT2(1),
    η=myT2(0),
    λ=myT2(1.113),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)
plot(x -> eigensol1d24.sol(x)[2], 1e-8, 0.002; label="λ1")

1/1.113

Optim.minimizer(opteigensol1)

eigensol1 = Eigensolve(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    rhomin=myT(0.00001),
    rho0=myT(10),
    A0=myT(Optim.minimizer(opteigensol1)[1]),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    λ=myT(Optim.minimizer(opteigensol1)[2]),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)

opteigensol2 = getlambda0(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    lambdalb=myT(-0.8),
    deltalambda=myT(0.2),
    rhomin=myT(0.00001),
    rho0=myT(10),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    rtol=1e-12,
    atol=1e-12,
    method=RadauIIA5(),
)

Optim.minimizer(opteigensol2)

eigensol2 = Eigensolve(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    rhomin=myT(0.00001),
    rho0=myT(10),
    A0=myT(Optim.minimizer(opteigensol2)[1]),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    λ=myT(Optim.minimizer(opteigensol2)[2]),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)

opteigensol3 = getlambda0(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    lambdalb=myT(-3.5),
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

Optim.minimizer(opteigensol3)

eigensol3 = Eigensolve(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    rhomin=myT(0.00001),
    rho0=myT(10),
    A0=myT(Optim.minimizer(opteigensol3)[1]),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    λ=myT(Optim.minimizer(opteigensol3)[2]),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)

opteigensol4 = getlambda0(
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

Optim.minimizer(opteigensol4)

eigensol4 = Eigensolve(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    rhomin=myT(0.00001),
    rho0=myT(10),
    A0=myT(Optim.minimizer(opteigensol4)[1]),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    λ=myT(Optim.minimizer(opteigensol4)[2]),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)

opteigensol5 = getlambda0(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    lambdalb=myT(-9.0),
    deltalambda=myT(1.0),
    rhomin=myT(0.00001),
    rho0=myT(10),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    rtol=1e-12,
    atol=1e-12,
    method=RadauIIA5(),
)

Optim.minimizer(opteigensol5)

eigensol5 = Eigensolve(
    Eigenfun;
    u1fun=itp1,
    u2fun=itp2,
    rhomin=myT(0.00001),
    rho0=myT(10),
    A0=myT(Optim.minimizer(opteigensol5)[1]),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    λ=myT(Optim.minimizer(opteigensol5)[2]),
    rtol=1e-16,
    atol=1e-16,
    method=RadauIIA5(),
)

plot(x -> eigensol1.sol(x)[1], 0.00001, 0.072; label="λ1")
plot!(x -> eigensol2.sol(x)[1], 0.00001, 0.072; label="λ2")
plot!(x -> eigensol3.sol(x)[1], 0.00001, 0.072; label="λ3")
plot!(x -> eigensol4.sol(x)[1], 0.00001, 0.072; label="λ4")
plot!(x -> eigensol5.sol(x)[1], 0.00001, 0.072; label="λ5")
eigensol1.sol(x)

writedlm("data/O1/d=3rho1data.dat", eigensol1.sol.t)
writedlm("data/O1/d=3u1data.dat", (collect(hcat(eigensol1.sol.u...)')))
writedlm("data/O1/d=3rho2data.dat", eigensol2.sol.t)
writedlm("data/O1/d=3u2data.dat", (collect(hcat(eigensol2.sol.u...)')))
writedlm("data/O1/d=3rho3data.dat", eigensol3.sol.t)
writedlm("data/O1/d=3u3data.dat", (collect(hcat(eigensol3.sol.u...)')))
writedlm("data/O1/d=3rho4data.dat", eigensol4.sol.t)
writedlm("data/O1/d=3u4data.dat", (collect(hcat(eigensol4.sol.u...)')))
writedlm("data/O1/d=3rho5data.dat", eigensol5.sol.t)
writedlm("data/O1/d=3u5data.dat", (collect(hcat(eigensol5.sol.u...)')))
writedlm(
    "data/O1/d=3parametersdata.dat",
    [Optim.minimizer(opteigensol1) Optim.minimizer(opteigensol2) Optim.minimizer(
        opteigensol3
    ) Optim.minimizer(opteigensol4) Optim.minimizer(opteigensol5)],
)
writedlm("data/O1/d=3fixsolt.dat", O1du0sol.t)
writedlm("data/O1/d=3fixsolu0.dat", O1du0sol.u1)
writedlm("data/O1/d=3fixsolu1.dat", O1du0sol.u2)
writedlm("data/O1/d=3fixsolu2.dat", O1du0sol.u3)
