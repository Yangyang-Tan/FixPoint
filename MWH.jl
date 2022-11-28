Udatasaver3!(
    usol_MWH;
    d=2.0,
    n=5.0,
    η=0.25,
    lb_A0=myT(0.006),
    ub_A0=0.046,
    lb_λ=0.01,
    ub_λ=1.8,
    rhomin_pre=0.1,
    tol_usol=1e-20,
    tol_A0=1e-16,
    tol_eigensol=1e-13,
    tol_λ=1e-13,
    dtmax=0.005,
)

Udatasaver3!(
    usol_MWH;
    d=2.0,
    n=1.0,
    η=0.15,
    lb_A0=myT(400.0),
    ub_A0=10^10,
    lb_λ=0.01,
    ub_λ=1.8,
    rhomin_pre=0.001,
    tol_usol=1e-20,
    tol_A0=1e-31,
    tol_eigensol=1e-13,
    tol_λ=1e-13,
    dtmax=0.005,
)

usol_MWH2[rationalize(0.15)]=usol_MWH[(2//1, 1//1)]

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
        usol_MWH2[rationalize(0.2)]|>etasave,
        usol_MWH2[rationalize(0.15)]|>etasave,
        usol_MWH2[rationalize(0.1)]|>etasave,
        usol_MWH2[rationalize(0.05)]|>etasave,
    ][:,4];
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
    [1.0, 0.5, 0.1, 0.05],
    [
        usol_lpap[2//1, rationalize(1.0)].η,
        usol_lpap[2//1, rationalize(1.5)].η,
        usol_lpap[2//1, rationalize(1.9)].η,
        usol_lpap[2//1, rationalize(1.95)].η,
    ];
    xaxis=:log,
    yaxis=:log,
    ylabel="η",
    xlabel="2-N",
    seriestype=:scatter,
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

usol_lpap[2//1, rationalize(1.9)].λ

usol_lpap[2//1, rationalize(1.95)].rhomin
