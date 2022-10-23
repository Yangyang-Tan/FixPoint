using JLD2
usol_WF = Dict()
# n=1
Udatasaver3!(
    usol_WF;
    d=2.1,
    n=1.0,
    lb_A0=myT(10),
    ub_A0=myT(1e13),
    lb_λ=0.3,
    ub_λ=1.5,
    rhomin_pre=1//100,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-16,
    tol_λ=1e-16,
)

Threads.@threads for i in 2.4:0.1:2.9
    Udatasaver3!(
        usol_WF;
        d=i,
        n=1.0,
        lb_A0=myT(10),
        ub_A0=myT(1e5),
        lb_λ=0.5,
        ub_λ=2.0,
        rhomin_pre=1//100,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end

Udatasaver3!(
    usol_WF;
    d=3.0,
    n=1.0,
    lb_A0=myT(1),
    ub_A0=myT(1000),
    lb_λ=0.1,
    ub_λ=2.0,
    rhomin_pre=1//100,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-16,
    tol_λ=1e-16,
)

Threads.@threads for i in 3.1:0.1:3.3
    Udatasaver3!(
        usol_WF;
        d=i,
        n=1.0,
        lb_A0=myT(1.0),
        ub_A0=myT(1000),
        lb_λ=0.5,
        ub_λ=2.0,
        rhomin_pre=1//100,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end

Threads.@threads for i in 3.4:0.1:3.9
    Udatasaver3!(
        usol_WF;
        d=i,
        n=1.0,
        lb_A0=myT(0.5),
        ub_A0=myT(1000),
        lb_λ=0.5,
        ub_λ=2.5,
        rhomin_pre=1//100,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end
save_object("data/UdataAll.jld2", usol_WF)
save_object("data/UdataAll_mini.jld2", usol_WF_mini)
minisave!(; udic_mini=usol_WF_mini, udic_large=usol_WF)
#n=2.0
Udatasaver3!(
    usol_WF;
    d=3.0,
    n=2.0,
    lb_A0=myT(0.5),
    ub_A0=myT(30),
    lb_λ=0.5,
    ub_λ=2.5,
    rhomin_pre=1//100,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-16,
    tol_λ=1e-16,
)

Threads.@threads for i in 3.4:0.1:3.9
    Udatasaver3!(
        usol_WF;
        d=i,
        n=2.0,
        lb_A0=myT(0.5),
        ub_A0=myT(1000),
        lb_λ=0.5,
        ub_λ=2.5,
        rhomin_pre=1//100,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end

#d=3.0

Udatasaver3!(
    usol_WF;
    d=3.0,
    n=10.0,
    lb_A0=myT(1.0),
    ub_A0=myT(3.0),
    lb_λ=0.5,
    ub_λ=2.0,
    rhomin_pre=1//100,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-16,
    tol_λ=1e-16,
)

Threads.@threads for i in 4.0:1.0:8.0
    Udatasaver3!(
        usol_WF;
        d=3.0,
        n=i,
        lb_A0=myT(1),
        ub_A0=myT(10),
        lb_λ=0.5,
        ub_λ=2.0,
        rhomin_pre=1//100,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end

Threads.@threads for i in 9.0:1.0:12.0
    Udatasaver3!(
        usol_WF;
        d=3.0,
        n=i,
        lb_A0=myT(0.2),
        ub_A0=myT(3.0),
        lb_λ=0.5,
        ub_λ=2.0,
        rhomin_pre=1//100,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end

#d=2.5

Udatasaver3!(
    usol_WF;
    d=2.5,
    n=10.0,
    lb_A0=myT2(0.187818),
    ub_A0=myT2(0.227818),
    lb_λ=0.2,
    ub_λ=1.5,
    rhomin_pre=1//100,
    tol_usol=1e-20,
    tol_A0=1e-18,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
)
Threads.@threads for i in 4.0:1.0:7.0
    Udatasaver3!(
        usol_WF;
        d=2.5,
        n=i,
        lb_A0=myT(0.5),
        ub_A0=myT(14.0),
        lb_λ=0.5,
        ub_λ=2.0,
        rhomin_pre=1//100,
        tol_usol=1e-21,
        tol_A0=1e-31,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end

#d=2.4

Udatasaver3!(
    usol_WF;
    d=2.4,
    n=10.0,
    lb_A0=myT(0.0428737),
    ub_A0=myT(0.0448647),
    lb_λ=0.2,
    ub_λ=0.8,
    rhomin_pre=1//100,
    tol_usol=1e-24,
    tol_A0=1e-30,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
)

#d=2.3

Udatasaver3!(
    usol_WF;
    d=2.3,
    n=10.0,
    lb_A0=myT(0.002),
    ub_A0=myT(0.0021),
    lb_λ=0.1,
    ub_λ=0.8,
    rhomin_pre=1//100,
    tol_usol=1e-20,
    tol_A0=1e-10,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
)

#d=2.2

Udatasaver3!(
    usol_WF;
    d=2.2,
    n=10.0,
    lb_A0=myT("1.41257442447300815170299538930733618e-06"),
    ub_A0=myT("1.41257442447308815170299538930733618e-06"),
    lb_λ=0.1,
    ub_λ=0.5,
    rhomin_pre=1//100,
    tol_usol=1e-24,
    tol_A0=1e-24,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
    dtmax=0.005,
)
#d=2.6
Udatasaver3!(
    usol_WF;
    d=2.6,
    n=3.0,
    lb_A0=myT(10),
    ub_A0=myT(50),
    lb_λ=0.2,
    ub_λ=2.0,
    rhomin_pre=1//10000,
    tol_usol=1e-20,
    tol_A0=1e-20,
    tol_eigensol=1e-16,
    tol_λ=1e-14,
    dtmax=0.005,
)

function testf(i)
    n = i - 2
    An = usol_WF[2.6, n].A0
    An1 = usol_WF[2.6, n + 1].A0
    mid =
        (2 + n)^((log(An) - log(An1)) / (log(n) - log(1 + n))) *
        exp((log(An1) * log(n) - log(An) * log(1 + n)) / (log(n) - log(1 + n)))
    lb = 0.75 * mid
    ub = 1.25 * mid
    Udatasaver3!(
        usol_WF;
        d=2.6,
        n=i,
        lb_A0=myT(lb),
        ub_A0=myT(ub),
        lb_λ=0.2,
        ub_λ=2.0,
        rhomin_pre=1//100,
        tol_usol=1e-24,
        tol_A0=1e-39,
        tol_eigensol=1e-16,
        tol_λ=1e-14,
        dtmax=0.005,
    )
    return (lb, ub)
end

testf(10.0)

#d=2.7
Udatasaver3!(
    usol_WF;
    d=2.7,
    n=3.0,
    lb_A0=myT(10),
    ub_A0=myT(100),
    lb_λ=0.2,
    ub_λ=2.0,
    rhomin_pre=1//10000,
    tol_usol=1e-20,
    tol_A0=1e-20,
    tol_eigensol=1e-16,
    tol_λ=1e-14,
    dtmax=0.005,
)

function WF_iter(d, i)
    n = i - 2
    An = usol_WF[2.6, n].A0
    An1 = usol_WF[2.6, n + 1].A0
    mid =
        (2 + n)^((log(An) - log(An1)) / (log(n) - log(1 + n))) *
        exp((log(An1) * log(n) - log(An) * log(1 + n)) / (log(n) - log(1 + n)))
    lb = 1.5 * mid
    ub = 4 * mid
    Udatasaver3!(
        usol_WF;
        d=d,
        n=i,
        lb_A0=myT(lb),
        ub_A0=myT(ub),
        lb_λ=0.2,
        ub_λ=1.5,
        rhomin_pre=1//1000,
        tol_usol=1e-24,
        tol_A0=1e-20,
        tol_eigensol=1e-16,
        tol_λ=1e-14,
        dtmax=0.005,
    )
    return (lb, ub)
end

for i in 10:10
    WF_iter(2.7, i)
end
#d=2.8

Threads.@threads for i in 9:1.0:10
    Udatasaver3!(
        usol_WF;
        d=2.8,
        n=i,
        lb_A0=myT(1),
        ub_A0=myT(4),
        lb_λ=0.5,
        ub_λ=1.5,
        rhomin_pre=1//10000,
        tol_usol=1e-24,
        tol_A0=1e-20,
        tol_eigensol=1e-16,
        tol_λ=1e-14,
        dtmax=0.005,
    )
end
Udatasaver3!(
    usol_WF;
    d=2.8,
    n=10.0,
    lb_A0=myT(1),
    ub_A0=myT(2),
    lb_λ=0.5,
    ub_λ=1.5,
    rhomin_pre=1//10000,
    tol_usol=1e-24,
    tol_A0=1e-20,
    tol_eigensol=1e-16,
    tol_λ=1e-14,
    dtmax=0.005,
)
#d=2.9
Udatasaver3!(
    usol_WF;
    d=2.9,
    n=7.0,
    lb_A0=myT(1),
    ub_A0=myT(15),
    lb_λ=0.5,
    ub_λ=1.5,
    rhomin_pre=1//10000,
    tol_usol=1e-20,
    tol_A0=1e-20,
    tol_eigensol=1e-16,
    tol_λ=1e-14,
    dtmax=0.005,
)
Threads.@threads for i in 4:1.0:7
    Udatasaver3!(
        usol_WF;
        d=2.9,
        n=i,
        lb_A0=myT(1),
        ub_A0=myT(15),
        lb_λ=0.5,
        ub_λ=1.5,
        rhomin_pre=1//10000,
        tol_usol=1e-20,
        tol_A0=1e-20,
        tol_eigensol=1e-16,
        tol_λ=1e-14,
        dtmax=0.005,
    )
end

#d=3.4
Udatasaver3!(
    usol_WF;
    d=3.4,
    n=5.0,
    lb_A0=myT(0.1),
    ub_A0=myT(5.0),
    lb_λ=0.5,
    ub_λ=1.999,
    rhomin_pre=1//10000,
    tol_usol=1e-20,
    tol_A0=1e-20,
    tol_eigensol=1e-16,
    tol_λ=1e-13,
    dtmax=0.005,
)


Threads.@threads for i in 3.4:0.1:3.9
    Udatasaver3!(
    usol_WF;
    d=i,
    n=10.0,
    lb_A0=myT(0.1),
    ub_A0=myT(5.0),
    lb_λ=0.5,
    ub_λ=1.999,
    rhomin_pre=1//10000,
    tol_usol=1e-20,
    tol_A0=1e-20,
    tol_eigensol=1e-16,
    tol_λ=1e-13,
    dtmax=0.005,
)
end
1
#n=1.5
Udatasaver3!(
    usol_WF;
    d=2.2,
    n=1.5,
    lb_A0=myT(usol_WF[2.2, 1.0].A0),
    ub_A0=myT(usol_WF[2.2, 2.0].A0),
    lb_λ=0.1,
    ub_λ=1.999,
    rhomin_pre=1//10000,
    tol_usol=1e-20,
    tol_A0=1e-20,
    tol_eigensol=1e-16,
    tol_λ=1e-13,
    dtmax=0.005,
)


Threads.@threads for i in 2.4:0.1:3.8
    Udatasaver3!(
    usol_WF;
    d=i,
    n=1.5,
    lb_A0=myT(usol_WF[i, 1.0].A0),
    ub_A0=myT(usol_WF[i, 2.0].A0),
    lb_λ=0.1,
    ub_λ=1.999,
    rhomin_pre=1//10000,
    tol_usol=1e-20,
    tol_A0=1e-20,
    tol_eigensol=1e-16,
    tol_λ=1e-13,
    dtmax=0.005,
)
end

usol_lpap[(3.0, 1.0)].λ

proc = run(`sleep 60`; wait=false)
Process(`sleep 60`, ProcessRunning)
kill(proc)
using Distributed


Threads.@threads for i in 1:8
    interrupt()
end

Udatasaver3!(
    usol_WF;
    d=3.2,
    n=4.0,
    lb_A0=myT(1.0),
    ub_A0=myT(7.0),
    lb_λ=0.5,
    ub_λ=1.9,
    rhomin_pre=1//10000,
    tol_usol=1e-20,
    tol_A0=1e-20,
    tol_eigensol=1e-16,
    tol_λ=1e-14,
    dtmax=0.005,
)
Threads.@threads for i in 4:1.0:9.0
    Udatasaver3!(
        usol_WF;
        d=i,
        n=10.0,
        lb_A0=myT(1.0),
        ub_A0=myT(10),
        lb_λ=0.5,
        ub_λ=1.9,
        rhomin_pre=1//10000,
        tol_usol=1e-20,
        tol_A0=1e-20,
        tol_eigensol=1e-16,
        tol_λ=1e-14,
        dtmax=0.005,
    )
end

1
for i in 10:10
    WF_iter(2.7, i)
end
#n=2.5
for i in 2.2:0.1:2.3
    Udatasaver2!(usol_WF; d=i, n=1.0, lb_A0=100.0, ub_A0=1e13, lb_λ=0.5, ub_λ=1.5)
end

tempv1 = 2.1:0.1:3.9
tempv2 = collect(2.1:0.1:3.9)

tempv12 = 2.2:0.1:3.9
tempv22 = collect(2.2:0.1:3.9)
tempv13 = 2.2:0.1:3.9
tempv23 = collect(2.2:0.1:3.9)

tempv1inf = 2.1:0.1:3.9
tempv2inf = 1 ./ (collect(2.1:0.1:3.9) .- 2)

tempv14 = 2.3:0.1:3.9
tempv24 = collect(2.3:0.1:4.9)
for i in eachindex(2.3:0.1:3.9)
    tempv24[i] = usol_WF[tempv14[i], 10.0].A0
end

for i in 2:6
    tempv24[i] = usol_WF[2.3, i].A0
end
tempv24[2:19]
plot(tempv24[2:6]; seriestype=:scatter, label="n=1", xaxis=:log10, yaxis=:log10)

for i in eachindex(2.1:0.1:3.9)
    tempv2[i] = usol_WF[tempv1[i], 1.0].λ
end

for i in eachindex(2.2:0.1:3.9)
    tempv22[i] = usol_WF[tempv12[i], 2.0].A0
end

plot(tempv1, tempv2; seriestype=:scatter, label="n=1")
plot(tempv12, tempv22; seriestype=:scatter, label="n=2")
plot(tempv13, tempv23; seriestype=:scatter, label="n=3")

plot(tempv1, 1 ./ tempv2; seriestype=:scatter, label="n=1")
plot!(tempv12, 1 ./ tempv22; seriestype=:scatter, label="n=2")
plot!(tempv13, 1 ./ tempv23; seriestype=:scatter, label="n=3")
plot!(tempv1inf, tempv2inf; seriestype=:scatter, label="n=inf")

# save_object("data/UdataAll.jld2", usol_WF)

usol_WF = load_object("data/UdataAll.jld2")
load("data/UdataAll.jld2", usol_WF[3.0, 2.0])

f = jldopen("data/UdataAll.jld2", "r")
close(ans)
f = jldopen("data/UdataAll.jld2", true, true, true, IOStream)

f["usol_WF"]

f

f["single_stored_object"]

@load "data/UdataAll.jld2" usol_WF["d", "n"]
tempsol = findUmin(;
    lb=myT(0.1),
    ub=myT(5e1),
    rhomin=myT(1e-2),
    rhomax=myT(10),
    d=myT(3.9),
    n=myT(2.0),
    η=myT(0.0),
    rtol=1e-20,
    atol=1e-20,
    method=Brent(),
)

for i in 3.6:0.1:3.9
    Udatasaver3!(
        usol_WF; d=i, n=10.0, lb_A0=0.1, ub_A0=10, lb_λ=1.0, ub_λ=2.0, rhomin_pre=1e-3
    )
end

for i in 3.0:0.1:3.5
    Udatasaver3!(
        usol_WF; d=i, n=3.0, lb_A0=0.1, ub_A0=40, lb_λ=0.5, ub_λ=2.0, rhomin_pre=1e-2
    )
end

for i in 2.9:0.1:2.9
    Udatasaver3!(
        usol_WF; d=i, n=3.0, lb_A0=1, ub_A0=30, lb_λ=0.2, ub_λ=2.0, rhomin_pre=1e-2
    )
end
usol_WF[2.9, 3.0]
for i in 3.3:0.1:3.7
    Udatasaver3!(
        usol_WF;
        d=i,
        n=10.0,
        lb_A0=myT(0.1),
        ub_A0=myT(2.1),
        lb_λ=1.0,
        ub_λ=1.95,
        rhomin_pre=1//100,
        tol_usol=1e-21,
        tol_A0=1e-30,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end

for i in 2.2:0.1:2.4
    Udatasaver3!(
        usol_WF; d=i, n=3.0, lb_A0=1, ub_A0=60, lb_λ=0.1, ub_λ=1.0, rhomin_pre=6e-1
    )
end
Udatasaver2!(
    usol_WF;
    d=3.0,
    n=10.0,
    lb_A0=myT2(
        "1.891036844326647060584417721890540596899559881481545401426986932414115402300479"
    ),
    ub_A0=myT2(
        "1.891036844326647060584417721890560596899559881481545401426986932414115402300479"
    ),
    lb_λ=1,
    ub_λ=1.5,
    rhomin_pre=1//10000000,
    tol_usol=1e-21,
    tol_A0=1e-40,
    tol_eigensol=1e-10,
    tol_λ=1e-5,
)

Udatasaver3!(
    usol_WF;
    d=2.3,
    n=8.0,
    lb_A0=myT(0.00986145859383776),
    ub_A0=myT(0.009898701698823532),
    lb_λ=0.1,
    ub_λ=0.8,
    rhomin_pre=1//100,
    tol_usol=1e-21,
    tol_A0=1e-30,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
)

Udatasaver3!(
    usol_WF;
    d=2.3,
    n=6.0,
    lb_A0=myT(0.08633006973034604),
    ub_A0=myT(0.08633006973234808),
    lb_λ=0.1,
    ub_λ=0.6,
    rhomin_pre=1//100,
    tol_usol=1e-21,
    tol_A0=1e-31,
    tol_eigensol=1e-14,
    tol_λ=1e-14,
)

Threads.@threads for i in 15.0:1.0:18.0
    Udatasaver3!(
        usol_WF;
        d=3.0,
        n=i,
        lb_A0=myT(0.5),
        ub_A0=myT(1.1),
        lb_λ=0.8,
        ub_λ=1.4,
        rhomin_pre=1//100,
        tol_usol=1e-21,
        tol_A0=1e-30,
        tol_eigensol=1e-16,
        tol_λ=1e-16,
    )
end

Udatasaver3!(
    usol_WF;
    d=3.0,
    n=8.0,
    lb_A0=myT(1.89),
    ub_A0=myT(5.0),
    lb_λ=0.9,
    ub_λ=1.4,
    rhomin_pre=1//100,
    tol_usol=1e-21,
    tol_A0=1e-30,
    tol_eigensol=1e-16,
    tol_λ=1e-16,
)

0.6351449110137829
0.5335461141427843
eltype([1.0, Double64(1.0), 1.0, 1, 1.01, 3, BigFloat(1.0)])
usol_WF[3.0, 10.0].A0

tempsol = du1solve(
    O1du1,
    myT(0.008611530822114632);
    rhomin=myT(1e-4),
    rhomax=myT(100),
    d=myT(3),
    n=myT(80.0),
    η=myT(0.0),
    rtol=1e-21,
    atol=1e-21,
    method=RadauIIA5(),
    dtmax=0.005,
)
1
plot(usol_WF[2.3, 1.0].t[1:10000], usol_WF[2.3, 1.0].u2[1:10000])

plot(tempsol.t[1:20000], tempsol.u3[1:20000])

tempeigensol = Eigensolve(
    Eigenfun;
    u1fun=Interpolations.interpolate((tempsol.t,), tempsol.u2, Gridded(Linear())),
    u2fun=Interpolations.interpolate((tempsol.t,), tempsol.u3, Gridded(Linear())),
    rhomin=tempsol.rhomin,
    rhomax=tempsol.rhomax,
    A0=myT(-1000000000000.0),
    d=tempsol.d,
    n=tempsol.n,
    η=tempsol.η,
    λ=myT(0.87),
    atol=1e-13,
    rtol=1e-13,
    method=RadauIIA5(),
)
tempsol.d
plot(tempeigensol.t[1:400], tempeigensol.u2[1:400])

dv1(;
    λ=0.5,
    v0=tempeigensol.u1[1],
    U1=Interpolations.interpolate((tempsol.t,), tempsol.u2, Gridded(Linear())),
    rhomin=1e-6,
    d=3,
    n=1,
    η=0,
)

tempsol.t
plot(tempsol.t[1:8000], tempsol.u1[1:8000])

tempsol.u1
tempsol.u1[1:18000]

plot(x -> x, -1, 1)
tempsol.u1[40000] - tempsol.u1[40001]
plot(tempsol.t[1:8000], tempsol.u2[1:8000])
