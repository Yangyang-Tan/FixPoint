usol_lpapsigma = Dict()

for i in 3.0:0.1:3.0
    d1p = rationalize(Float32(i))
    n1p = 1//1
    Udatasaver3!(
        usol_lpapsigma;
        d=usol_lpap[d1p, n1p].d,
        n=usol_lpap[d1p, n1p].n,
        η=etafunsigma(usol_lpap[d1p, n1p], i) * 0.16,
        lb_A0=usol_lpap[d1p, n1p].A0 * 0.1,
        ub_A0=usol_lpap[d1p, n1p].A0 * 1.5,
        lb_λ=0.1,
        ub_λ=1.99,
        rhomin_pre=usol_lpap[d1p, n1p].rhomin,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-14,
        tol_λ=1e-14,
        dtmax=0.005,
    )
end

for i in 2.0:0.1:2.0
    d1p = rationalize(Float32(i))
    n1p = 1//1
    Udatasaver3!(
        usol_lpapsigma;
        d=usol_lpapsigma[d1p, n1p].d,
        n=usol_lpapsigma[d1p, n1p].n,
        η=(
            0.2 * etafunsigma(usol_lpapsigma[d1p, n1p], i) +
            0.8 * usol_lpapsigma[d1p, n1p].η
        ),
        lb_A0=usol_lpapsigma[d1p, n1p].A0 * 0.5,
        ub_A0=usol_lpapsigma[d1p, n1p].A0 * 10.5,
        lb_λ=0.8,
        ub_λ=2.5,
        rhomin_pre=usol_lpapsigma[d1p, n1p].rhomin,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-14,
        tol_λ=1e-14,
        dtmax=0.005,
    )
end

etafunsigma(usol_lpap[3//1, 1//1], 2//1)
etafunsigma(usol_lpapsigma[2//1, 1//1], 2//1)
usol_lpapsigma[2//1, 1//1].η

etafunsigma(usol_lpapsigma[3//1, 1//1], 3//1)
etasave("data/lpap/d=3.0_n=1_sigma.dat", usol_lpapsigma[3//1, 1//1])
1

plot([i for i in 2.1:0.1:3.9], [usol_lpap[i//10, 4//1].λ for i in 21:39])
writedlm("data/lpap/nu_lpa_prime_n=4_y.dat", 1 ./ [usol_lpap[i//10, 4//1].λ for i in 21:39])
writedlm("data/lpap/nu_lpa_prime_n=4_x.dat", [i for i in 2.1:0.1:3.9])

writedlm("data/lpap/nu_lpa_n=4_y.dat", 1 ./ [usol_WF[i / 10, 4.0].λ for i in 22:39])
writedlm("data/lpap/nu_lpa_n=4_x.dat", [i for i in 2.2:0.1:3.9])
