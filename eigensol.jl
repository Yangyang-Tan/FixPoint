eigensol1 = Eigensolve(
    Eigenfun,
    u1fun=itp1,
    u2fun=itp2,
    rhomin=myT(0.00001),
    rho0=myT(10),
    A0=myT(2.0),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    λ=myT(1.62),
    rtol=1e-20,
    atol=1e-20,
    method=RadauIIA5(),
)
plot(x -> eigensol1.sol(x)[2], 0.00001, 0.07, label="λ1")

eigensol2 = Eigensolve(
    Eigenfun,
    u1fun=itp1,
    u2fun=itp2,
    rhomin=myT(0.00001),
    rho0=myT(10),
    A0=myT(242000.0),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    λ=myT(val_lambda2),
    rtol=1e-20,
    atol=1e-20,
    method=RadauIIA5(),
)

eigensol5 = Eigensolve(
    Eigenfun,
    u1fun=itp1,
    u2fun=itp2,
    rhomin=myT(0.00001),
    rho0=myT(10),
    A0=myT(-415200000000000.0),
    d=myT(3.0),
    n=myT(1),
    η=myT(0),
    λ=myT(val_lambda5),
    rtol=1e-20,
    atol=1e-20,
    method=RadauIIA5(),
)


plot(x -> eigensol1.sol(x)[2], 0.00001, 0.07, label="λ1")


plot(x -> eigensol5.sol(x)[1], 0.00001, 0.00005, label="λ1")

plot(x -> eigensol1.sol(x)[1] - eigensol1.sol(0.00001)[1], 0.00001, 0.07, label="λ1")
plot!(x -> eigensol2.sol(x)[1] - eigensol2.sol(0.00001)[1], 0.00001, 0.07, label="λ2")
plot!(x -> eigensol3.sol(x)[1] - eigensol3.sol(0.00001)[1], 0.00001, 0.07, label="λ3")
plot!(x -> eigensol4.sol(x)[1] - eigensol4.sol(0.00001)[1], 0.00001, 0.07, label="λ4")
plot!(x -> eigensol5.sol(x)[1] - eigensol5.sol(0.00001)[1], 0.00001, 0.08, label="λ5")


plot(x -> eigensol1.sol(x)[2] - eigensol1.sol(0.00001)[2], 0.00001, 0.07, label="λ1")
plot!(x -> eigensol2.sol(x)[2] - eigensol2.sol(0.00001)[2], 0.00001, 0.07, label="λ2")
plot!(x -> eigensol3.sol(x)[2] - eigensol3.sol(0.00001)[2], 0.00001, 0.07, label="λ3")
plot!(x -> eigensol4.sol(x)[2] - eigensol4.sol(0.00001)[2], 0.00001, 0.07, label="λ4")
plot!(x -> eigensol5.sol(x)[2] - eigensol5.sol(0.00001)[2], 0.00001, 0.07, label="λ5")




1

plot!(x -> eigensol4.sol(x)[1] - eigensol4.sol(0.00001)[1], 0.00001, 0.8, label="λ4")
plot!(x -> eigensol5.sol(x)[1] - eigensol5.sol(0.00001)[1], 0.00001, 0.8, label="λ5")
