plot(O1du0sol.t[1:1000], O1du0sol.u2[1:1000])
plot(O1sol.t[1:6000], O1sol.u3[1:6000])

plot(O1sol.t[1:10], O1sol.u3[1:10]; label="u''(ρ)", title="d=2.4,N=1")
plot!(
    O1solerr1.t[1:2000], O1solerr1.u3[1:2000]; label="u''(ρ),A=A0+1e-12", title="d=2.4,N=1"
)
plot!(
    O1solerr2.t[1:2000], O1solerr2.u3[1:2000]; label="u''(ρ),A=A0+1e-12", title="d=2.4,N=1"
)

Nderivative(O1sol.t[1:5000], O1sol.u3[1:5000], 2 * 1e-5)

tempsol = du1solve(
    O1du1,
    myT(usol_WF[2.2,8.0].A0);
    rhomin=myT(1//10),
    rhomax=myT(100),
    d=myT(2.2),
    n=myT(8.0),
    η=myT(0.0),
    rtol=1e-26,
    atol=1e-26,
    method=RadauIIA5(),
    dtmax=0.005,
)
usol_WF[2.2,8.0]
plot(t -> tempsol.sol(t)[2], 1e-1, 2.0)



plot(t -> usol_WF[2.2,8.0].sol(t)[2], 1e-1, 2.0)

usol_WF[2.2,10.0].rhomin
usol_WF[2.2, 10.0].λ
usol_WF[3.0, 3.0].A0
usol_WF[3.0, 8.0].A0
showall(usol_WF[2.2, 10.0].A0)
usol_WF[2.5, 4.0].A0
0.2046251118182223
myT("4.35025400087700902945657626460808982e-04")
usol_WF[2.2, 9.0].λ
0.6351449370833719
[1/usol_WF[d, n].λ for d in 2.2:0.1:3.9, n in 1.0:1.0:10.0]

[0.6351449370833719, 0.2046251118182223]
