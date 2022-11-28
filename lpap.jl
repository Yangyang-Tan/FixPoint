using JLD2
usol_lpap = Dict()

usol_MWH2 = Dict()
usol_WF = load_object("data/UdataAll.jld2")
save_object("data/UlpadataAll.jld2", usol_lpap)

usol_lpap = load_object("data/UlpadataAll.jld2")
usol_lpap2

for i in 1:24
    usol_lpap[rationalize.(Float32.(collect(keys(usol_lpap2))[i]))] = usol_lpap2[collect(
        keys(usol_lpap2)
    )[i]]
end

usol_lpap

usol_lpap2.keys[8] === (2.1, 1.0)
usol_lpap2.keys[15] == (2.3, 1.0)

usol_lpap[2.3, 1.0] = usol_lpap[usol_lpap.keys[15]]

etafun(rho0, upp, d)
etafun(usol_WF[2.2, 2.0], 2.2)

usol_WF[3.0, 1.0].d
usol_WF[2.2, 1.0].A0
usol_WF[2.5, 1.0].rhomax
Udatasaver3!(
    usol_lpap;
    d=usol_WF[2.5, 1.0].d,
    n=usol_WF[2.5, 1.0].n,
    η=etafun(usol_WF[2.5, 1.0], 2.5),
    lb_A0=usol_WF[2.5, 1.0].A0 * 0.01,
    ub_A0=usol_WF[2.5, 1.0].A0 * 1.5,
    lb_λ=0.2,
    ub_λ=1.99,
    rhomin_pre=1//10000,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    dtmax=0.005,
)

etafun(usol_lpap[2.5, 1.0], 2.5)
etafun(usol_lpap[2.1, 1.0], 2.1)
etafun(usol_lpap[2.2, 1.0], 2.2)

Udatasaver3!(
    usol_lpap;
    d=usol_lpap[d1p, n1p].d,
    n=usol_lpap[d1p, n1p].n,
    η=etafun(usol_lpap[d1p, n1p], i),
    lb_A0=usol_lpap[d1p, n1p].A0 * 0.1,
    ub_A0=usol_lpap[d1p, n1p].A0 * 10.5,
    lb_λ=0.95,
    ub_λ=1.99,
    rhomin_pre=usol_lpap[d1p, n1p].rhomin,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    dtmax=0.005,
)

Threads.@threads for i in 2.2:0.1:2.2
    d1p = rationalize(Float32(i))
    n1p = 1//1
    Udatasaver3!(
        usol_lpap;
        d=usol_lpap[d1p, n1p].d,
        n=usol_lpap[d1p, n1p].n,
        η=etafun(usol_lpap[d1p, n1p], i),
        lb_A0=usol_lpap[d1p, n1p].A0 * 0.1,
        ub_A0=usol_lpap[d1p, n1p].A0 * 10.5,
        lb_λ=0.95,
        ub_λ=1.99,
        rhomin_pre=usol_lpap[d1p, n1p].rhomin,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-14,
        tol_λ=1e-14,
        dtmax=0.005,
    )
end

Threads.@threads for i in 3.0:0.1:3.9
    d1p = rationalize(Float32(i))
    n1p = 1//1
    Udatasaver3!(
        usol_lpap;
        d=usol_WF[i, 2.0].d,
        n=usol_WF[i, 2.0].n,
        η=etafun(usol_WF[i, 2.0], i),
        lb_A0=usol_WF[i, 2.0].A0 * 0.5,
        ub_A0=usol_WF[i, 2.0].A0 * 3,
        lb_λ=0.2,
        ub_λ=1.99,
        rhomin_pre=usol_WF[i, 2.0].rhomin,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-14,
        tol_λ=1e-14,
        dtmax=0.005,
    )
end
@qthreads for i in 3.2:0.1:3.5
    d1p = rationalize(Float32(i))
    n1p = 2//1
    Udatasaver3!(
        usol_lpap;
        d=usol_lpap[d1p, n1p].d,
        n=usol_lpap[d1p, n1p].n,
        η=etafun(usol_lpap[d1p, n1p], i),
        lb_A0=usol_lpap[d1p, n1p].A0 * 0.5,
        ub_A0=usol_lpap[d1p, n1p].A0 * 2,
        lb_λ=usol_WF[max(i, 2.2), 2.0].λ * 0.5,
        ub_λ=2.0,
        rhomin_pre=usol_lpap[d1p, n1p].rhomin,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-14,
        tol_λ=1e-14,
        dtmax=0.005,
    )
end

usol_lpap[3//1, 1//1].η

Threads.@threads for x in 1:20
    println("$x $(Threads.threadid())")
end

Threads.@threads for i in 3.1:0.1:3.9
    d1p = rationalize(Float32(i))
    n1p = 1//1
    Udatasaver3!(
        usol_lpap;
        d=usol_lpap[d1p, n1p].d,
        n=usol_lpap[d1p, n1p].n,
        η=etafun(usol_lpap[d1p, n1p], i),
        lb_A0=usol_lpap[d1p, n1p].A0 * 0.1,
        ub_A0=usol_lpap[d1p, n1p].A0 * 10.5,
        lb_λ=1.1,
        ub_λ=2.5,
        rhomin_pre=usol_lpap[d1p, n1p].rhomin,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
        dtmax=0.005,
    )
end

for i in 1:5
    Threads.@threads for i in 2.0:0.1:2.0
        d1p = rationalize(Float32(i))
        n1p = 1//1
        Udatasaver3!(
            usol_lpap;
            d=usol_lpap[d1p, n1p].d,
            n=usol_lpap[d1p, n1p].n,
            η=etafun(usol_lpap[d1p, n1p], i),
            lb_A0=usol_lpap[d1p, n1p].A0 * 0.1,
            ub_A0=usol_lpap[d1p, n1p].A0 * 10.5,
            lb_λ=0.5,
            ub_λ=1.99,
            rhomin_pre=usol_lpap[d1p, n1p].rhomin,
            tol_usol=1e-21,
            tol_A0=1e-31,
            tol_eigensol=1e-14,
            tol_λ=1e-14,
            dtmax=0.005,
        )
    end
end

Threads.@threads for i in 2.0:0.1:2.0
    d1p = rationalize(Float32(i))
    n1p = 1//1
    Udatasaver3!(
        usol_lpap;
        d=usol_lpap[d1p, n1p].d,
        n=usol_lpap[d1p, n1p].n,
        η=etafun(usol_lpap[d1p, n1p], i),
        lb_A0=usol_lpap[d1p, n1p].A0 * 0.9,
        ub_A0=usol_lpap[d1p, n1p].A0 * 1.1,
        lb_λ=0.5,
        ub_λ=1.8,
        rhomin_pre=usol_lpap[d1p, n1p].rhomin,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
        dtmax=0.005,
    )
end

vlpap = [1 / usol_lpap[rationalize(i), 1//1].λ for i in 2.0:0.1:3.9]
vlpa = [1 / usol_WF[i, 1.0].λ for i in 2.1:0.1:3.9]

A0vlpap = [usol_lpap[rationalize(i), 1//1].A0 for i in 2.0:0.1:3.9]

veta = [usol_lpap[rationalize(i), 1//1].η for i in 2.0:0.1:3.9]

writedlm("data/vlpap_n=1.dat", vlpap)

writedlm("data/vlpa_n=1.dat", vlpa)
writedlm("data/veta_n=1.dat", veta)

plot(2.1:0.1:3.9, vlpa; seriestype=:scatter)
plot!(2.1:0.1:3.9, vlpap; seriestype=:scatter)
plot(0.1:0.1:1.9, A0vlpap; seriestype=:scatter, xaxis=:log, yaxis=:log)

Udatasaver3!(
    usol_lpap;
    d=usol_lpap[3.9, 1.0].d,
    n=usol_lpap[3.9, 1.0].n,
    η=etafun(usol_lpap[3.9, 1.0], 3.9),
    lb_A0=usol_lpap[3.9, 1.0].A0 * 0.1,
    ub_A0=usol_lpap[3.9, 1.0].A0 * 100.5,
    lb_λ=0.3,
    ub_λ=1.99,
    rhomin_pre=usol_lpap[3.9, 1.0].rhomin,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    dtmax=0.005,
)

for i in 1:5
    Udatasaver3!(
        usol_lpap;
        d=usol_lpap[2.1, 1.0].d,
        n=usol_lpap[2.1, 1.0].n,
        η=etafun(usol_lpap[2.1, 1.0], 2.1),
        lb_A0=usol_lpap[2.1, 1.0].A0 * 0.5,
        ub_A0=usol_lpap[2.1, 1.0].A0 * 2.5,
        lb_λ=0.3,
        ub_λ=1.7,
        rhomin_pre=1//10000,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
        dtmax=0.005,
    )
end

testv = rand(19, 6)

etasave(usol_lpap[39//10, 1//1])

for i in 21:1:39
    testv[i - 20, :] .= etasave(usol_lpap[i//10, 4//1])
end
testv
usol_lpap[20//10, 1//1].d

writedlm("data/lpap/n=4.dat", testv)

tempsol = usol_lpap[3.0, 1.0].sol
tempsol = du1solve(
    O1du1,
    myT(usol_TF[2.9999, 1.0].A0);
    rhomin=myT(1//10000),
    rhomax=myT(10),
    d=myT(2.9999),
    n=myT(1.0),
    η=myT(0.0),
    rtol=1e-20,
    atol=1e-20,
    method=RadauIIA5(),
    dtmax=0.005,
)
plot(usol_WF[3.0, 1.0].t, usol_WF[3.0, 1.0].u1)
plot(usol_lpap[22//10, 4//1].t[1:5000], usol_lpap[22//10, 4//1].u2[1:5000])
plot(usol_WF[2.1, 1.0].t[1:10000], usol_WF[2.1, 1.0].u2[1:10000])
usol_WF[2.1, 1.0].λ
usol_lpap[2.2, 1.0].λ
usol_lpap

usol_lpap
