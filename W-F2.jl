usol_WF2 = Dict()
Udatasaver3!(
    usol_WF2;
    d=3.9,
    n=1.0,
    lb_A0=myT(0.0000001),
    ub_A0=myT(25.23392888283688),
    lb_λ=1.0,
    ub_λ=1.999,
    rhomin_pre=1//1000000,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-16,
    tol_λ=1e-16,
)

Udatasaver3!(
    usol_WF2;
    d=3.999999,
    n=1.0,
    lb_A0=myT(0.0000001),
    ub_A0=myT(5.23392888283688),
    lb_λ=1.0,
    ub_λ=1.9999999,
    rhomin_pre=1//10000000,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-16,
    tol_λ=1e-16,
)

plot(x->usol_WF2[(3.9, 1.0)].sol(x)[1], 1e-5,2e-2)
plot(x -> usol_WF2[(3.99, 1.0)].sol(x)[1], 1e-5, 2e-2)
plot(x -> usol_WF2[(3.999, 1.0)].sol(x)[1], 1e-5, 2e-2)
plot(x -> usol_WF2[(3.999999, 1.0)].sol(x)[1], 1e-5, 2e-2)
8/3
save_object("data/UdataAll_WF2.jld2", usol_WF2)
