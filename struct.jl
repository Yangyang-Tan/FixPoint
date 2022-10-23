struct USolution{T}
    sol::Any
    t::AbstractArray
    u1::AbstractArray
    u2::AbstractArray
    u3::AbstractArray
    A0::T
    d::T
    n::T
    η::T
    rhomin::T
    rhomax::T
    λ
end

struct USolution_mini{T}
    A0::T
    d::T
    n::T
    η::T
    rhomin::T
    rhomax::T
    λ
end

function USolution(sol, t, u1, u2)
    return USolution(
        sol,
        t,
        u1,
        u2,
        similar(t),
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

function USolution(sol, t, u1, u2, u3)
    return USolution(
        sol, t, u1, u2, u3, nothing, nothing, nothing, nothing, nothing, nothing, nothing
    )
end

function USolution(sol, t, u1, u2, A0, d, n, η, rhomin, rhomax, λ)
    return USolution(sol, t, u1, u2, similar(t), A0, d, n, η, rhomin, rhomax, λ)
end

function USolution(
    sol::T;
    A0=nothing,
    d=nothing,
    n=nothing,
    η=nothing,
    rhomin=nothing,
    rhomax=nothing,
    λ=nothing,
) where {T<:OrdinaryDiffEq.AbstractDiffEqArray}
    # @cast v2[i][j] := sol.u[j][i]
    v2 = reverse.(invert(sol.u))
    return USolution(sol, reverse(sol.t), v2..., A0, d, n, η, rhomin, rhomax, λ)
end

function USolution_mini(u::USolution)
    return USolution_mini(
        u.A0,
        u.d,
        u.n,
        u.η,
        u.rhomin,
        u.rhomax,
        u.λ,
    )
end
