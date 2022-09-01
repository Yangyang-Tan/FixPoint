#Flow eqn for O(1) and u'(ρ)
function O1du1(du, u, p, rho; d::T, n::T, η::T = 0) where {T}
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    du[1] = u[2]
    du[2] = u[3]
    du[3] =
        -(3 * u[2] + 2 * rho * u[3]) / (6 * pi^2 * (1 + u[1] + 2 * rho * u[2])^2) +
        rho * u[2] - 2 * u[1]
end

#Flow eqn for O(1) and u(ρ), u[1] for u(ρ), u[2] for u'(ρ), u[3] for u''(ρ)
function O1du0(du, u, p, rho)
    # c1 =(2+d-η)/(pi^(d/2)*(d*(2+d)*gamma(d/2)*2^(d-1)))
    d, n, η, c1 = p
    du[1] = u[2]
    du[2] = u[3]
    du[3] =
        c1 * ((-1 + n) / (1 + u[2]) + (1 + u[2] + 2 * rho * u[3])^(-1)) - (d * u[1]) +
        (-2 + d + η) * rho * u[2]
end

#Eigen perturbation ode at O(ϵ)
function Eigenfun(du, u, p, rho, U1, U2)
    d, n, η, λ = p
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    # @show rho
    du[1] = u[2]
    du[2] =
        ((-d + λ) * u[1] * (1 + U1(rho) + 2 * rho * U2(rho))^2) / (2 * c1 * rho) +
        (1 / (2 * c1 * rho)) *
        u[2] *
        (
            rho * (-2 + d + η) * (1 + U1(rho) + 2 * rho * U2(rho))^2 -
            (c1 * (-1 + n) * (1 + U1(rho) + 2 * rho * U2(rho))^2) / (1 + U1(rho))^2 - c1
        )
end



function test1()
    A0 = 809.0
    rho = 10.0
    rho0 = rho
    u = [
        1 / (150 * A0 * pi^2 * rho0^2) + (A0 * rho0^3) / 3,
        A0 * rho0^2 - (1 / (75 * A0 * pi^2 * rho0^3)),
        2 * A0 * rho0 + 1 / (25 * A0 * pi^2 * rho0^4),
    ]
    # du[1] = u[2]
    # du[2] = u[3]
    # du[3] = 1 / (6 * pi^2 * (1 + u[2] + 2 * rho * u[3])) + rho * u[2] - 3 * u[1]
    return 1 / (6 * pi^2 * (1 + u[2] + 2 * rho * u[3])) + rho * u[2] - 3 * u[1]
end

function geterror(sol)
    @. -(3 * sol.u2 + 2 * sol.t * sol.u3) /
       (6 * pi^2 * (1 + sol.u1 + 2 * sol.t * sol.u2)^2) + sol.t * sol.u2 - 2 * sol.u1
end

function dv1(λ, U1, rho, d = 3, n = 1, η = 0)
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    return (c1^-1) * (λ - d) * (1 + U1(rho))^2 / n
end
