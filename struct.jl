struct USolution
    sol::Any
    t::AbstractArray
    u1::AbstractArray
    u2::AbstractArray
    u3::AbstractArray
end

function USolution(sol::T) where {T <: AbstractDiffEqArray}
    # @cast v2[i][j] := sol.u[j][i]
    v2 = reverse.(invert(sol.u))
    USolution(sol, reverse(sol.t), v2...)
end

function USolution(sol, t, u1, u2)
    USolution(sol, t, u1, u2, similar(t))
end
