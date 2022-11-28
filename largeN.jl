for i in 3.0:0.1:3.0
    Udatasaver3!_LN(
        usol_WF;
        d=i,
        n=40.0,
        lb_A0=myT(0.12802119650749802),
        ub_A0=myT(0.13002119650749802),
        lb_λ=0.5,
        ub_λ=2.5,
        rhomin_pre=1//50,
        tol_usol=1e-28,
        tol_A0=1e-31,
        tol_eigensol=1e-14,
        tol_λ=1e-14,
    )
end
