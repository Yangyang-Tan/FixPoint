usol_lpap[(22//10, 4//1)].A0
plot(22:1:29, [usol_lpap[(i//10, 2//1)].A0 for i in 22:1:29])
plot(usol_lpap[22//10, 4//1].t[1:5000], usol_lpap[22//10, 4//1].u2[1:5000])
etafun(usol_lpap[22//10, 4//1], 22//10)
plot!([etafun(usol_lpap[i//10, 4//1], i//10) for i in 22:1:29])
plot!([etafun(usol_WF[i / 10, 4//1], i / 10) for i in 22:1:29])

usol_lpap[(23//10, 3//1)].A0 - usol_lpap[(22//10, 3//1)].A0

plot(usol_lpap[21//10, 3//1].t[1:10000], usol_lpap[21//10, 3//1].u2[1:10000])
usol_lpap[21//10, 4//1].rhomin
Udatasaver3!(
    usol_lpap;
    d=usol_lpap[21//10, 2//1].d,
    n=usol_lpap[21//10, 2//1].n,
    η=etafun(usol_lpap[21//10, 2//1], 21//10),
    lb_A0=usol_lpap[21//10, 2//1].A0 * 0.5,
    ub_A0=usol_lpap[21//10, 2//1].A0 * 2.1,
    lb_λ=0.05,
    ub_λ=0.8,
    rhomin_pre=1 / 100,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    dtmax=0.005,
)

Udatasaver3!(
    usol_lpap;
    d=usol_lpap[21//10, 2//1].d,
    n=usol_lpap[21//10, 2//1].n,
    η=etafun(usol_lpap[21//10, 2//1], 21//10),
    lb_A0=usol_lpap[21//10, 2//1].A0 * 0.5,
    ub_A0=usol_lpap[21//10, 2//1].A0 * 2.1,
    lb_λ=0.05,
    ub_λ=0.8,
    rhomin_pre=1 / 100,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    dtmax=0.005,
)

Udatasaver3!(
    usol_lpap;
    d=41//20,
    n=usol_lpap[21//10, 4//1].n,
    η=usol_lpap[21//10, 4//1].η,
    lb_A0=usol_lpap[21//10, 4//1].A0 * 0.5,
    ub_A0=usol_lpap[21//10, 4//1].A0 * 2.1,
    lb_λ=0.01,
    ub_λ=0.2,
    rhomin_pre=1 / 10,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    dtmax=0.005,
)

etasave(usol_WF[3.0, 4.0])
etasave(usol_lpap[30//10, 4//1])
1 / usol_lpap[20//10, 1//1].λ

plot(usol_lpap[22//10, 4//1].t[1:10000], usol_lpap[22//10, 4//1].u2[1:10000])

plot!(21:1:39, [usol_lpap[(i//10, 2//1)].η for i in 21:1:39])

plot(21:1:39, [1 / usol_lpap[(i//10, 4//1)].λ for i in 21:1:39])

plot(usol_lpapsigma[20//10, 1//1].t[1:5000], usol_lpapsigma[20//10, 1//1].u2[1:5000])
plot(usol_WF[30//10, 100//1].t[1:5000], usol_WF[30//10, 100//1].u2[1:5000])

plot!(x -> usol_WF[30//10, 30//1].sol(x)[2], 0.1, 1)

plot!(x -> usol_WF[30//10, 40//1].sol(x)[2], 0.1, 1)
1

usol_WF[39//10, 100//1]
usol_WF[30//10, 40//1].A0
usol_WF[3.0, 9.0].A0

1 / usol_WF[30//10, 40//1].λ

1 / 1.019595273624346

plot([etasave(usol_lpap[i//10, 4//1])[4] for i in 21:1:39])

plot([etasave(usol_lpap[i//10, 4//1])[4] for i in 21:1:39])
[1 / usol_WF[i, 10.0].λ for i in 2.2:0.1:3.9]
writedlm(
    "/home/tyy/Documents/FixPoint/data/lpap/nu_lpa_prime_n=10_y.dat",
    [1 / usol_WF[i, 10.0].λ for i in 2.2:0.1:3.9],
)

plot(x -> usol_lpap[20//10, 19//10].sol(x)[2], 0.0001, 1.2)
plot(x -> usol_lpap[20//10, rationalize(1.95)].sol(x)[2], 0.1, 2.1)

plot(x -> usol_CS[20//10, rationalize(1.0)].sol(x)[2], 0.1, 1.4)
etafun_CS(usol_CS[20//10, rationalize(1.0)], 20//10)
usol_CS[20//10, rationalize(1.0)].A0

plot(x -> usol_CS[20//10, rationalize(0.8)].sol(x)[2], 0.01, 1.1)
etafun_CS(usol_CS[20//10, rationalize(0.8)], 20//10)
usol_CS[20//10, rationalize(0.8)].A0
usol_CS[20//10, rationalize(0.9)]
testint = interpolate(myT.(randn(10)), BSpline(Linear()))
cuitp = adapt(Vector{Double64}, testint);

myT(1.0):myT(1.0):myT(10.0)
myT.(1.0:1.0:10.0)
cuitp(myT.(1.0:1.0:10.0))
unsafe_trunc(Int64, x::Double64)=trunc(Int64, x)
