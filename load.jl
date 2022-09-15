using DifferentialEquations, Plots, Roots, DoubleFloats, ModelingToolkit
using SpecialFunctions, Dierckx, SplitApplyCombine,BSplineKit, Interpolations,DelimitedFiles
using Optim

function Nderivative(x::AbstractArray, y::AbstractArray, z::Number, derivative::Int=2, order::Int=6)
    spl = interpolate(x, y, BSplineOrder(order))
    D2f = diff(spl, Derivative(derivative))
    return D2f(z)
end
