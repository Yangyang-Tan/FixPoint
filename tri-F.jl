usol_TF = Dict()
1
@time using Plots
@time using OrdinaryDiffEq
1
using Plots
using CUDA
Udatasaver3!(
    usol_TF;
    d=2.999999,
    n=1.0,
    lb_A0=myT(0.0000001),
    ub_A0=myT(25.23392888283688),
    lb_λ=0.3,
    ub_λ=1.5,
    rhomin_pre=1//10000,
    tol_usol=1e-20,
    tol_A0=1e-31,
    tol_eigensol=1e-16,
    tol_λ=1e-16,
)
save_object("data/UdataAll_T2.jld2", usol_TF)
findeigenval2_y2(
    Eigenfun;
    Usol=usol_WF[2.2, 9.0],
    lambdalb=-2.5,
    lambdaub=-0.000001,
    tol_eigensol=1e-14,
    tol_λ=1e-29,
    method=Brent(),
)
plot(usol_WF[2.2, 8.0].t[1:50000], usol_WF[2.2, 8.0].u2[1:50000])

eigensol_WF = Dict()
eigensol_T2 = Dict()
eigensol_T3 = Dict()
eigensol_T4 = Dict()

eigensol_WF[2.4, 1.0, 1] = findeigenval3(
    Eigenfun;
    Usol=usol_WF[2.4, 1.0],
    lambdalb=0.1,
    lambdaub=2.0,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)
eigensol_WF[2.4, 1.0, 2] = findeigenval3(
    Eigenfun;
    Usol=usol_WF[2.4, 1.0],
    lambdalb=-1.0,
    lambdaub=-0.1,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)
eigensol_WF[2.4, 1.0, 3] = findeigenval3(
    Eigenfun;
    Usol=usol_WF[2.4, 1.0],
    lambdalb=-4.0,
    lambdaub=-1.0,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)
eigensol_WF[2.4, 1.0, 4] = findeigenval3(
    Eigenfun;
    Usol=usol_WF[2.4, 1.0],
    lambdalb=-7.0,
    lambdaub=-4.0,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)
eigensol_WF[2.4, 1.0, 5] = findeigenval3(
    Eigenfun;
    Usol=usol_WF[2.4, 1.0],
    lambdalb=-9.0,
    lambdaub=-6.5,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)
plot(x -> x, -1, 1)
usol_WF[2.4, 1.0].rhomin
plot(
    x -> eigensol_WF[2.4, 1.0, 1].sol(x)[1] / eigensol_WF[2.4, 1.0, 1].sol(1e-3)[1],
    1e-3,
    0.135,
)
plot!(
    x -> eigensol_WF[2.4, 1.0, 2].sol(x)[1] / eigensol_WF[2.4, 1.0, 2].sol(1e-3)[1],
    1e-3,
    0.135,
)
plot!(
    x -> eigensol_WF[2.4, 1.0, 3].sol(x)[1] / eigensol_WF[2.4, 1.0, 3].sol(1e-3)[1],
    1e-3,
    0.135,
)
plot!(
    x -> eigensol_WF[2.4, 1.0, 4].sol(x)[1] / eigensol_WF[2.4, 1.0, 4].sol(1e-3)[1],
    1e-3,
    0.135,
)
plot!(
    x -> eigensol_WF[2.4, 1.0, 5].sol(x)[1] / eigensol_WF[2.4, 1.0, 5].sol(1e-3)[1],
    1e-3,
    0.135,
)

usol_T2 = Dict()
usol_T3 = Dict()
usol_T4 = Dict()

Udatasaver3!(
    usol_T2;
    d=2.4,
    n=1.0,
    lb_A0=myT(
        "1.043463300189046733047927848991565776584122946731511177716190679303755426633109"
    ),
    ub_A0=myT(
        "1.043463500189046733047927848991565776584122946731511177716190679303755426633109"
    ),
    lb_λ=0.003,
    ub_λ=1.5,
    rhomin_pre=1//1000000,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
)
eigensol_T2[2.4, 1.0, 1] = findeigenval3(
    Eigenfun;
    Usol=usol_T2[2.4, 1.0],
    lambdalb=1.8,
    lambdaub=2.3,
    tol_eigensol=1e-21,
    tol_λ=1e-20,
    method=RadauIIA5(),
)

eigensol_T2[2.4, 1.0, 2] = findeigenval3(
    Eigenfun;
    Usol=usol_T2[2.4, 1.0],
    lambdalb=0.1,
    lambdaub=1.0,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)

eigensol_T2[2.4, 1.0, 3] = findeigenval3(
    Eigenfun;
    Usol=usol_T2[2.4, 1.0],
    lambdalb=-1.0,
    lambdaub=-0.1,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)

eigensol_T2[2.4, 1.0, 4] = findeigenval3(
    Eigenfun;
    Usol=usol_T2[2.4, 1.0],
    lambdalb=-2.2,
    lambdaub=-1.0,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)

eigensol_T2[2.4, 1.0, 5] = findeigenval3(
    Eigenfun;
    Usol=usol_T2[2.4, 1.0],
    lambdalb=-4.2,
    lambdaub=-2.2,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)

plot(
    x -> eigensol_T2[2.4, 1.0, 1].sol(x)[1] / eigensol_T2[2.4, 1.0, 1].sol(1e-3)[1],
    1e-3,
    0.5,
)
plot!(
    x -> eigensol_T2[2.4, 1.0, 2].sol(x)[1] / eigensol_T2[2.4, 1.0, 2].sol(1e-3)[1],
    1e-3,
    0.5,
)
plot!(
    x -> eigensol_T2[2.4, 1.0, 3].sol(x)[1] / eigensol_T2[2.4, 1.0, 3].sol(1e-3)[1],
    1e-3,
    0.5,
)
plot!(
    x -> eigensol_T2[2.4, 1.0, 4].sol(x)[1] / eigensol_T2[2.4, 1.0, 4].sol(1e-3)[1],
    1e-3,
    0.5,
)
plot!(
    x -> eigensol_T2[2.4, 1.0, 5].sol(x)[1] / eigensol_T2[2.4, 1.0, 5].sol(1e-3)[1],
    1e-3,
    0.5,
)

Udatasaver3!(
    usol_T3;
    d=2.4,
    n=1.0,
    lb_A0=myT(
        "0.006650286507215782001088138419950468947830290362523216311050394728775002213255677",
    ),
    ub_A0=myT(
        "0.006650286527235784001088138419950468947830290362523216311050394728775002213255677",
    ),
    lb_λ=0.003,
    ub_λ=1.5,
    rhomin_pre=1//1000000,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
)

eigensol_T3[2.4, 1.0, 1] = findeigenval3(
    Eigenfun;
    Usol=usol_T3[2.4, 1.0],
    lambdalb=1.8,
    lambdaub=2.1,
    tol_eigensol=1e-21,
    tol_λ=1e-20,
    method=RadauIIA5(),
)
eigensol_T3[2.4, 1.0, 2] = findeigenval3(
    Eigenfun;
    Usol=usol_T3[2.4, 1.0],
    lambdalb=1.0,
    lambdaub=1.7,
    tol_eigensol=1e-17,
    tol_λ=1e-17,
    method=RadauIIA5(),
)

eigensol_T3[2.4, 1.0, 3] = findeigenval3(
    Eigenfun;
    Usol=usol_T3[2.4, 1.0],
    lambdalb=0.1,
    lambdaub=1.0,
    tol_eigensol=1e-17,
    tol_λ=1e-15,
    method=RadauIIA5(),
)

eigensol_T3[2.4, 1.0, 4] = findeigenval3(
    Eigenfun;
    Usol=usol_T3[2.4, 1.0],
    lambdalb=-0.1,
    lambdaub=-1.0,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)

eigensol_T3[2.4, 1.0, 5] = findeigenval3(
    Eigenfun;
    Usol=usol_T3[2.4, 1.0],
    lambdalb=-1.5,
    lambdaub=-1.0,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    method=RadauIIA5(),
)

plot(
    x -> eigensol_T3[2.4, 1.0, 1].sol(x)[1] / eigensol_T3[2.4, 1.0, 1].sol(1e-3)[1],
    1e-3,
    1.1,
)
plot!(
    x -> eigensol_T3[2.4, 1.0, 2].sol(x)[1] / eigensol_T3[2.4, 1.0, 2].sol(1e-3)[1],
    1e-3,
    1.1,
)
plot!(
    x -> eigensol_T3[2.4, 1.0, 3].sol(x)[1] / eigensol_T3[2.4, 1.0, 3].sol(1e-3)[1],
    1e-3,
    1.1,
)
plot!(
    x -> eigensol_T3[2.4, 1.0, 4].sol(x)[1] / eigensol_T3[2.4, 1.0, 4].sol(1e-3)[1],
    1e-3,
    1.1,
)
plot!(
    x -> eigensol_T3[2.4, 1.0, 5].sol(x)[1] / eigensol_T3[2.4, 1.0, 5].sol(1e-3)[1],
    1e-3,
    1.1,
)

Udatasaver3!(
    usol_T4;
    d=2.4,
    n=1.0,
    lb_A0=myT(
        "0.0001012417040816608903347432443585211143212138860398738948839017969019409495128915",
    ),
    ub_A0=myT(
        "0.0001012457040816608903347432443585211143212138860398738948839017969019409495128915",
    ),
    lb_λ=0.003,
    ub_λ=1.5,
    rhomin_pre=1//1000000,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
)

eigensol_T2[2.4, 1.0, 1].λ
eigensol_T3[2.4, 1.0, 1].λ

eigensol_T4[2.4, 1.0, 1] = findeigenval3(
    Eigenfun;
    Usol=usol_T4[2.4, 1.0],
    lambdalb=1.9,
    lambdaub=2.3,
    tol_eigensol=1e-21,
    tol_λ=1e-23,
    method=RadauIIA5(),
)
eigensol_T4[2.4, 1.0, 2] = findeigenval3(
    Eigenfun;
    Usol=usol_T4[2.4, 1.0],
    lambdalb=1.4,
    lambdaub=1.7,
    tol_eigensol=1e-21,
    tol_λ=1e-17,
    method=RadauIIA5(),
)

eigensol_T4[2.4, 1.0, 3] = findeigenval3(
    Eigenfun;
    Usol=usol_T4[2.4, 1.0],
    lambdalb=0.8,
    lambdaub=1.2,
    tol_eigensol=1e-21,
    tol_λ=1e-16,
    method=RadauIIA5(),
)

eigensol_T4[2.4, 1.0, 4] = findeigenval3(
    Eigenfun;
    Usol=usol_T4[2.4, 1.0],
    lambdalb=0.1,
    lambdaub=0.8,
    tol_eigensol=1e-21,
    tol_λ=1e-16,
    method=RadauIIA5(),
)

eigensol_T4[2.4, 1.0, 5] = findeigenval3(
    Eigenfun;
    Usol=usol_T4[2.4, 1.0],
    lambdalb=-0.5,
    lambdaub=-0.01,
    tol_eigensol=1e-15,
    tol_λ=1e-14,
    method=RadauIIA5(),
)

plot(
    x -> eigensol_T4[2.4, 1.0, 1].sol(x)[1] / eigensol_T4[2.4, 1.0, 1].sol(1e-3)[1],
    1e-3,
    1.8,
)
plot!(
    x -> eigensol_T4[2.4, 1.0, 2].sol(x)[1] / eigensol_T4[2.4, 1.0, 2].sol(1e-3)[1],
    1e-3,
    1.8,
)
plot!(
    x -> eigensol_T4[2.4, 1.0, 3].sol(x)[1] / eigensol_T4[2.4, 1.0, 3].sol(1e-3)[1],
    1e-3,
    1.8,
)
plot!(
    x -> eigensol_T4[2.4, 1.0, 4].sol(x)[1] / eigensol_T4[2.4, 1.0, 4].sol(1e-3)[1],
    1e-3,
    1.8,
)
plot!(
    x -> eigensol_T4[2.4, 1.0, 5].sol(x)[1] / eigensol_T4[2.4, 1.0, 5].sol(1e-3)[1],
    1e-3,
    1.8,
)
save_object("data/O1/d=2.4/eigensol_WF.jld2", eigensol_WF)
save_object("data/O1/d=2.4/eigensol_T2.jld2", eigensol_T2)
save_object("data/O1/d=2.4/eigensol_T3.jld2", eigensol_T3)
save_object("data/O1/d=2.4/eigensol_T4.jld2", eigensol_T4)
