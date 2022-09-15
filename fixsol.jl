val_A0
@time O1du0sol2 = du0solve(
    O1du0, myT(val2_A0); rtol=1e-18, atol=1e-18, method=RadauIIA5(), d=myT(2.5))

val2_A0



plot(O1du0sol2.t[1:200], O1du0sol2.u3[1:200])

val2_A0 = find_zero(
    x -> getA0(x; d=myT(2.5)),
    (myT(val2_A0-0.001), myT(val2_A0+0.001)),
    Bisection();
    rtol=1e-4,
    atol=1e-4,
)




writedlm("data/O1/d=2.5fixsol2t.dat", O1du0sol2.t)
writedlm("data/O1/d=2.5fixsol2u0.dat", O1du0sol2.u1)
writedlm("data/O1/d=2.5fixsol2u1.dat", O1du0sol2.u2)
writedlm("data/O1/d=2.5fixsol2u2.dat", O1du0sol2.u3)




writedlm("data/O1/d=2.5fixsolt.dat", O1du0sol2.t)
writedlm("data/O1/d=2.5fixsolu0.dat", O1du0sol2.u1)
writedlm("data/O1/d=2.5fixsolu1.dat", O1du0sol2.u2)
writedlm("data/O1/d=2.5fixsolu2.dat", O1du0sol2.u3)



showall(val_A0)
plot(O1du0sol2.t[1:1000], O1du0sol2.u1[1:1000])

@time O1du0sol3 = du0solve(
    O1du0, myT2(0.0003); rtol=1e-12, atol=1e-12, method=RadauIIA5(), d=myT2(2.1)
)
