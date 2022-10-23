
writedlm("data/A0/A0.dat", [usol_WF[d, n].A0 for d in 2.2:0.1:3.9, n in 1.0:1.0:10.0])

writedlm("data/A0/A0_n=1.5.dat", [usol_WF[d, 1.5].A0 for d in 2.2:0.1:3.9])

writedlm("data/nu/nu.dat", [1 / usol_WF[d, n].λ for d in 2.2:0.1:3.9, n in 1.0:1.0:10.0])

writedlm("data/nu/nu_n=1.5.dat", [1 / usol_WF[d, 1.5].λ for d in 2.2:0.1:3.9])


usol_WF[2.3, 10.0].A0

plot(usol_WF[2.2, 8.0].t[1:33000], usol_WF[2.2, 8.0].u2[1:33000])

plot([usol_WF[2.3, i].λ for i in 2:8])
[usol_WF[3.0, i].λ for i in 2:8]
usol_WF[2.3, 6].
