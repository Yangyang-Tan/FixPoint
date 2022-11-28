using OrdinaryDiffEq, Plots, Roots, DoubleFloats
using SpecialFunctions, Dierckx, SplitApplyCombine, Interpolations, DelimitedFiles
using Optim, ThreadPools

# function Nderivative(
#     x::AbstractArray, y::AbstractArray, z::Number, derivative::Int=2, order::Int=6
# )
#     spl = interpolate(x, y, BSplineOrder(order))
#     D2f = diff(spl, Derivative(derivative))
#     return D2f(z)
# end
