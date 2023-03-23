using ThreadPools
#usol_MWH is LPA' and η is computed iteratively
#########################################################
ndata_Opt =
    rationalize.([
        1.95,
        1.94,
        1.93,
        1.92,
        1.91,
        1.89,
        1.88,
        1.86,
        1.83,
        1.81,
        1.78,
        1.74,
        1.7,
        1.65,
        1.59,
        1.53,
        1.45,
        1.36,
        1.26,
        1.14,
        1.0,
    ])

etadata_Opt =
    rationalize.([
        0.03,
        0.033,
        0.036,
        0.04,
        0.044,
        0.048,
        0.053,
        0.058,
        0.064,
        0.07,
        0.077,
        0.085,
        0.094,
        0.103,
        0.113,
        0.124,
        0.137,
        0.15,
        0.165,
        0.182,
        0.2,
    ])
#ini
let d = rationalize(2.0), n = rationalize(2.0), η = 0.1
    Udatasaver3!(
        usol_MWH_eta;
        d=d,
        n=n,
        η=η,
        lb_A0= 0.9*myT(0.43073977664059515),
        ub_A0=2.0000000000000001 * 0.43073977664059515,
        lb_λ=0.0001,
        ub_λ=0.2,
        rhomin_pre=0.001,
        tol_usol=1e-20,
        tol_A0=1e-31,
        tol_eigensol=1e-12,
        tol_λ=1e-12,
        dtmax=0.005,
    )
    et = rationalize(η)
    usol_MWH_d2n2_eta[et] = usol_MWH_eta[d, n]
end



let d = rationalize(2.0), n = rationalize(2.0), η = 0.03
    Udatasaver3!(
        usol_MWH_eta;
        d=d,
        n=n,
        η=η,
        lb_A0=0.99999999999999 * myT(2.918380700552567e-31),
        ub_A0=1.0000000000000001 * 2.918380700552567e-31,
        lb_λ=0.002,
        ub_λ=0.4,
        rhomin_pre=0.1,
        tol_usol=1e-28,
        tol_A0=1e-31,
        tol_eigensol=1e-23,
        tol_λ=1e-23,
        dtmax=0.005,
    )
    et = rationalize(η)
    usol_MWH_d2n2_eta[et] = usol_MWH_eta[d, n]
end

#iter
for i in 1:1
    Udatasaver3!(
        usol_MWH;
        d=2.0,
        n=usol_MWH[2//1, rationalize(1.94)].n,
        η=0.0 * usol_MWH[2//1, rationalize(1.94)].η +
          1 * etafun(usol_MWH[2//1, rationalize(1.94)], 2.0),
        lb_A0=0.8 * usol_MWH[2//1, rationalize(1.94)].A0,
        ub_A0=1.5 * usol_MWH[2//1, rationalize(1.94)].A0,
        lb_λ=0.001,
        ub_λ=0.1,
        rhomin_pre=0.01,
        tol_usol=1e-20,
        tol_A0=1e-31,
        tol_eigensol=1e-12,
        tol_λ=1e-12,
        dtmax=0.005,
    )
end
etafun(usol_MWH[2//1, rationalize(1.94)], 2.0)
usol_MWH[2//1, rationalize(1.94)].η
#qbthreads
#@qbthreads
Threads.@threads for j in ndata_Opt[2:7]
    for i in 1:5
        Udatasaver3!(
            usol_MWH;
            d=2.0,
            n=usol_MWH[2//1, j].n,
            η=0.0 * usol_MWH[2//1, j].η + 1.0 * etafun(usol_MWH[2//1, j], 2.0),
            lb_A0=0.7 * usol_MWH[2//1, j].A0,
            ub_A0=2.0 * usol_MWH[2//1, j].A0,
            lb_λ=0.01,
            ub_λ=0.5,
            rhomin_pre=0.01,
            tol_usol=1e-18,
            tol_A0=1e-31,
            tol_eigensol=1e-10,
            tol_λ=1e-10,
            dtmax=0.005,
        )
    end
end

Threads.@threads for j in ndata_Opt[2:8]
    for i in 1:1
        Udatasaver3!(
            usol_MWH;
            d=2.0,
            n=usol_MWH[2//1, j].n,
            η=0.0 * usol_MWH[2//1, j].η + 1.0 * etafun(usol_MWH[2//1, j], 2.0),
            lb_A0=0.8 * usol_MWH[2//1, j].A0,
            ub_A0=1.5 * usol_MWH[2//1, j].A0,
            lb_λ=0.001,
            ub_λ=0.1,
            rhomin_pre=0.01,
            tol_usol=1e-27,
            tol_A0=1e-31,
            tol_eigensol=1e-22,
            tol_λ=1e-22,
            dtmax=0.005,
        )
    end
end

GC.gc(true)

ndata_Opt
plot(
    2 .- ndata_Opt,
    (x -> 1/etasave(usol_MWH[2//1, x])[4]).(ndata_Opt);
    xaxis=:log,
    yaxis=:log,
    seriestype=:scatter,
    ylabel="1/ρ_0",
    xlabel="2-N"
)



etafun(usol_MWH[2//1, 195//100], 2.0)

etafun(usol_MWH[2//1, 195//100], 2.0)
1
etafun(usol_MWH[2//1, 195//100], 2.0)
Udatasaver3!(
    usol_MWH;
    d=2.0,
    n=usol_MWH[2//1, 195//100].n,
    η=0.2 * usol_MWH[2//1, 195//100].η + 0.8 * etafun(usol_MWH[2//1, 195//100], 2.0),
    lb_A0=0.1 * usol_MWH[2//1, 195//100].A0,
    ub_A0=1.2 * usol_MWH[2//1, 195//100].A0,
    lb_λ=0.002,
    ub_λ=0.5,
    rhomin_pre=0.2,
    tol_usol=1e-18,
    tol_A0=1e-31,
    tol_eigensol=1e-12,
    tol_λ=1e-12,
    dtmax=0.005,
)

for k in 1:6
    Threads.@threads for j in [150//100, 190//100, 195//100]
        Udatasaver3!(
            usol_MWH;
            d=2.0,
            n=usol_MWH[2//1, j].n,
            η=0.4 * usol_MWH[2//1, j].η + 0.6 * etafun(usol_MWH[2//1, j], 2.0),
            lb_A0=0.4 * usol_MWH[2//1, j].A0,
            ub_A0=2.0 * usol_MWH[2//1, j].A0,
            lb_λ=0.002,
            ub_λ=0.5,
            rhomin_pre=0.1,
            tol_usol=1e-18,
            tol_A0=1e-31,
            tol_eigensol=1e-12,
            tol_λ=1e-12,
            dtmax=0.005,
        )
    end
end

for k in 1:4
    Threads.@threads for j in [150//100, 190//100, 195//100]
        Udatasaver3!(
            usol_MWH;
            d=2.0,
            n=usol_MWH[2//1, j].n,
            η=0.2 * usol_MWH[2//1, j].η + 0.8 * etafun(usol_MWH[2//1, j], 2.0),
            lb_A0=0.1 * usol_MWH[2//1, j].A0,
            ub_A0=2.0 * usol_MWH[2//1, j].A0,
            lb_λ=0.002,
            ub_λ=0.5,
            rhomin_pre=0.1,
            tol_usol=1e-18,
            tol_A0=1e-31,
            tol_eigensol=1e-12,
            tol_λ=1e-12,
            dtmax=0.005,
        )
    end
end

for k in 1:6
    Threads.@threads for j in [150//100, 190//100, 195//100]
        Udatasaver3!(
            usol_MWH;
            d=2.0,
            n=usol_MWH[2//1, j].n,
            η=0.2 * usol_MWH[2//1, j].η + 0.8 * etafun(usol_MWH[2//1, j], 2.0),
            lb_A0=0.1 * usol_MWH[2//1, j].A0,
            ub_A0=2.0 * usol_MWH[2//1, j].A0,
            lb_λ=0.002,
            ub_λ=0.5,
            rhomin_pre=0.1,
            tol_usol=1e-18,
            tol_A0=1e-31,
            tol_eigensol=1e-12,
            tol_λ=1e-12,
            dtmax=0.005,
        )
    end
end

etafun(usol_MWH[2//1, 198//100], 2.0)
1 / usol_MWH[(2//1, 199//100)].λ
usol_MWH2[rationalize(0.15)] = usol_MWH[(2//1, 1//1)]

plot(
    1 ./ [4.0, 5.0, 6.0, 7.0],
    [
        usol_MWH[2//1, rationalize(4.0)].λ,
        usol_MWH[2//1, rationalize(5.0)].λ,
        usol_MWH[2//1, rationalize(6.0)].λ,
        usol_MWH[2//1, rationalize(7.0)].λ,
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="1/ν",
    xlabel="1/N",
    label="",
    seriestype=:scatter,
)

plot(
    [0.2, 0.15, 0.1, 0.05],
    [
        etasave(usol_MWH2[rationalize(0.2)]),
        etasave(usol_MWH2[rationalize(0.15)]),
        etasave(usol_MWH2[rationalize(0.1)]),
        etasave(usol_MWH2[rationalize(0.05)]),
    ][
        :, 4
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="1/ν",
    xlabel="η",
    label="",
    seriestype=:scatter,
)

plot!(
    [0.2, 0.15, 0.1, 0.05],
    [
        usol_MWH2[rationalize(0.2)].λ,
        usol_MWH2[rationalize(0.15)].λ,
        usol_MWH2[rationalize(0.1)].λ,
        usol_MWH2[rationalize(0.05)].λ,
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="1/ν",
    xlabel="η",
    label="",
    # seriestype=:scatter,
)

plot!(
    1 ./ [4.0, 5.0, 6.0, 7.0],
    [
        usol_MWH[2//1, rationalize(4.0)].λ,
        usol_MWH[2//1, rationalize(5.0)].λ,
        usol_MWH[2//1, rationalize(6.0)].λ,
        usol_MWH[2//1, rationalize(7.0)].λ,
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="1/ν",
    xlabel="1/N",
    label="",
    # seriestype=:scatter,
)

plot(
    1 ./ [4.0, 5.0, 6.0, 7.0],
    1 ./ [
        etasave(usol_MWH[2//1, rationalize(4.0)])[4],
        etasave(usol_MWH[2//1, rationalize(5.0)])[4],
        etasave(usol_MWH[2//1, rationalize(6.0)])[4],
        etasave(usol_MWH[2//1, rationalize(7.0)])[4],
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="1/ρ_0",
    xlabel="1/N",
    label="",
    seriestype=:scatter,
)

usol_lpap[20//10, 1//1].A0

for j in 1:4
    for i in 2.0:0.1:2.0
        d1p = rationalize(Float32(i))
        n1p = rationalize(1.95)
        Udatasaver3!(
            usol_lpap;
            d=usol_lpap[d1p, n1p].d,
            n=usol_lpap[d1p, n1p].n,
            η=etafun(usol_lpap[d1p, n1p], i),
            lb_A0=usol_lpap[d1p, n1p].A0 * 0.3,
            ub_A0=usol_lpap[d1p, n1p].A0 * 1.5,
            lb_λ=0.001,
            ub_λ=0.1,
            rhomin_pre=usol_lpap[d1p, n1p].rhomin,
            tol_usol=1e-21,
            tol_A0=1e-16,
            tol_eigensol=1e-13,
            tol_λ=1e-13,
            dtmax=0.005,
        )
    end
end

usol_lpap[2//1, 2//2].rhomin
etafun(usol_lpap[2//1, rationalize(1.95)], 2)

plot(
    [1.0, 0.5, 0.1, 0.05, 0.02],
    [
        usol_MWH[2//1, rationalize(1.0)].η,
        usol_MWH[2//1, rationalize(1.5)].η,
        usol_MWH[2//1, rationalize(1.9)].η,
        usol_MWH[2//1, rationalize(1.95)].η,
        usol_MWH[2//1, rationalize(1.98)].η,
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="η",
    xlabel="2-N",
    seriestype=:scatter,
)

plot!(
    [1.0, 0.5, 0.1, 0.05, 0.02],
    [
        usol_MWH[2//1, rationalize(1.0)].η,
        usol_MWH[2//1, rationalize(1.5)].η,
        usol_MWH[2//1, rationalize(1.9)].η,
        usol_MWH[2//1, rationalize(1.95)].η,
        etafun(usol_MWH[2//1, 198//100], 2.0),
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="η",
    xlabel="2-N",
    # seriestype=:scatter,
)

plot!(
    [1.0, 0.5, 0.1, 0.05],
    [
        usol_lpap[2//1, rationalize(1.0)].λ,
        usol_lpap[2//1, rationalize(1.5)].λ,
        usol_lpap[2//1, rationalize(1.9)].λ,
        usol_lpap[2//1, rationalize(1.95)].λ,
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="1/ν",
    xlabel="2-N",
    label="",
    # seriestype=:scatter,
)

plot(
    [1.0, 0.5, 0.1, 0.05],
    [
        usol_lpap[2//1, rationalize(1.0)].λ,
        usol_lpap[2//1, rationalize(1.5)].λ,
        usol_lpap[2//1, rationalize(1.9)].λ,
        usol_lpap[2//1, rationalize(1.95)].λ,
    ];
    xaxis=:log,
    yaxis=:log,
    seriestype=:scatter,
)

etasave(usol)

plot(
    [1.0, 0.5, 0.1, 0.05],
    [
        etasave(usol_lpap[2//1, rationalize(1.0)])[4],
        etasave(usol_lpap[2//1, rationalize(1.5)])[4],
        etasave(usol_lpap[2//1, rationalize(1.9)])[4],
        etasave(usol_lpap[2//1, rationalize(1.95)])[4],
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="ρ_0",
    xlabel="2-N",
    seriestype=:scatter,
)

plot(
    [1.0, 0.5, 0.1, 0.05],
    [
        usol_lpap[2//1, rationalize(1.0)].A0,
        usol_lpap[2//1, rationalize(1.5)].A0,
        usol_lpap[2//1, rationalize(1.9)].A0,
        usol_lpap[2//1, rationalize(1.95)].A0,
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="1/ν",
    xlabel="2-N",
    label="",
    # seriestype=:scatter,
)

USolution_mini(usol_MWH[2//1, rationalize(1.94)])
for i in ndata_Opt
    usol_MWH_mini[2//1, i] = USolution_mini(usol_MWH[2//1, i])
end

usol_MWH_mini
save_object("data/MWH/usol_mini/usol_MWH_Opt_mini.jld2", usol_MWH_mini)
save_object("data/MWH/usol_full/usol_MWH_Opt_full.jld2", usol_MWH)
