using DifferentialEquations, Optimization, OptimizationPolyalgorithms, OptimizationOptimJL,SciMLSensitivity, Zygote, Plots
function lossEigen(p)
    sol = solve(prob, Tsit5(), p=p, saveat=tsteps)
    loss = sum(abs2, sol .- 1)
    return loss, sol
end

1
