#Flow eqn for O(1) and u'(ρ)
function O1du1(du, u, p, rho)
    d, n, η, c1 = p
    # c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    du[1] = u[2]
    du[2] = u[3]
    du[3] =
        (
            (1 + u[2] + 2 * rho * u[3])^2 * ((η - 2) * u[2] + (η - 2 + d) * rho * u[3]) -
            c1 * u[3] * (3 + (n - 1) * (1 + (2 * rho * u[3]) / (1 + u[2]))^2)
        ) / (2 * c1 * rho)
    return nothing
end
#Flow eqn for C-S Regulator and u'(ρ)
#u[1] for u'(ρ), u[2] for u''(ρ), du[2] for u'''(ρ)
function effective_CS(du, u, p, rho)
    d, n, η, c1 = p
    # c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    c3 = (2 - η) * gamma(2 - d / 2) / (2 * (4 * π)^(d / 2))
    du[1] = u[2]
    du[2] = u[3]
    du[3] =
        -1 / 2 * (
            -((-2 + η) * u[2] * (1 + u[2] + 2 * rho * u[3])^(2 - d / 2)) +
            u[3] * (
                3 * c3 - rho * (-2 + d + η) * (1 + u[2] + 2 * rho * u[3])^(2 - d / 2) +
                c3 * (-1 + n) * (1 + (2 * rho * u[3]) / (1 + u[2]))^(2 - d / 2)
            )
        ) / (c3 * rho)
    # @show du[2] =(-3*c2*u(2)+(-2*u(1)+rho*u(2))*(1+u(1)+2*rho*u(2))^2)/(2*c2*rho)
    return nothing
end

function effective_rfun(d, nr, η)
    lp = lpfun(d, nr, η)
    return function effective_r(du, u, p, rho)
        d, n, η, c1 = p
        # c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
        c4 = (1 / 2) * (1 / (4 * pi)^(d / 2)) * (1 / gamma(d / 2))
        du[1] = u[2]
        du[2] = u[3]
        du[3] =
            (-3 * u[3]) / (2 * rho) -
            ((-2 + η) * u[2] + (rho * (-2 + d + η) + c4 * (-1 + n) * lp(u[2])) * u[3]) /
            (2 * c4 * rho * lp(u[2] + 2 * rho * u[3]))
        return nothing
    end
end

# function effective_r(du, u, p, rho)
#     d, n, η, c1 = p
#     # c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
#     c4 = (1 / 2) * (1 / (4 * pi)^(d / 2)) * (1 / gamma(d / 2))
#     du[1] = u[2]
#     du[2] = u[3]
#     du[3] =
#         (-3 * u[3]) / (2 * rho) -
#         ((-2 + η) * u[2] + (rho * (-2 + d + η) + c4 * (-1 + n) * lp(u[2])) * u[3]) /
#         (2 * c4 * rho * lp(u[2] + 2 * rho * u[3]))
#     # @show du[2] =(-3*c2*u(2)+(-2*u(1)+rho*u(2))*(1+u(1)+2*rho*u(2))^2)/(2*c2*rho)
#     return nothing
# end

function effective_CS(du, u, p, rho)
    d, n, η, c1 = p
    # c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    c3 = (2 - η) * gamma(2 - d / 2) / (2 * (4 * π)^(d / 2))
    du[1] = u[2]
    du[2] = u[3]
    du[3] =
        -1 / 2 * (
            -(
                (-2 + η) *
                u[2] *
                sign(1 + u[2] + 2 * rho * u[3]) *
                (abs(1 + u[2] + 2 * rho * u[3]))^(2 - d / 2)
            ) +
            u[3] * (
                3 * c3 -
                rho *
                (-2 + d + η) *
                sign(1 + u[2] + 2 * rho * u[3]) *
                (abs(1 + u[2] + 2 * rho * u[3]))^(2 - d / 2) +
                c3 *
                (-1 + n) *
                sign(1 + (2 * rho * u[3]) / (1 + u[2])) *
                (abs(1 + (2 * rho * u[3]) / (1 + u[2])))^(2 - d / 2)
            )
        ) / (c3 * rho)
    # @show du[2] =(-3*c2*u(2)+(-2*u(1)+rho*u(2))*(1+u(1)+2*rho*u(2))^2)/(2*c2*rho)
    return nothing
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
    return nothing
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
    return nothing
end

#Eigen perturbation ode at O(ϵ) for C-S Regulator
#u[1] for v(ρ), u[2] for v'(ρ), du[2] for v''(ρ)
function Eigenfun_CS(du, u, p, rho, U1, U2)
    d, n, η, λ = p
    c3 = (2 - η) * gamma(2 - d / 2) / (2 * (4 * π)^(d / 2))
    # @show rho
    du[1] = u[2]
    du[2] =
        (
            -(c3 * u[2]) +
            ((-d + λ) * u[1] + rho * (-2 + d + η) * u[2]) *
            (abs(1 + U1(rho) + 2 * rho * U2(rho)))^(2 - d / 2) -
            c3 *
            (-1 + n) *
            u[2] *
            (abs(1 + (2 * rho * U2(rho)) / (1 + U1(rho))))^(2 - d / 2)
        ) / (2 * c3 * rho)
    return nothing
end

# function Eigenfun_r(du, u, p, rho, U1, U2)
#     d, n, η, λ = p
#     c4 = (1 / 2) * (1 / (4 * pi)^(d / 2)) * (1 / gamma(d / 2))
#     # @show rho
#     du[1] = u[2]
#     du[2] =
#         (
#             (d - λ) * u[1] -
#             (
#                 rho * (-2 + d + η) +
#                 c4 * (-1 + n) * lp(U1(rho)) +
#                 c4 * lp(U1(rho) + 2 * rho * U2(rho))
#             ) * u[2]
#         ) / (2 * c4 * rho * lp(U1(rho) + 2 * rho * U2(rho)))
#     return nothing
# end

function Eigenfun_rfun(d, nr, η)
    lp = lpfun(d, nr, η)
    return function Eigenfun_r(du, u, p, rho, U1, U2)
        d, n, η, λ = p
        c4 = (1 / 2) * (1 / (4 * pi)^(d / 2)) * (1 / gamma(d / 2))
        # @show rho
        du[1] = u[2]
        du[2] =
            (
                (d - λ) * u[1] -
                (
                    rho * (-2 + d + η) +
                    c4 * (-1 + n) * lp(U1(rho)) +
                    c4 * lp(U1(rho) + 2 * rho * U2(rho))
                ) * u[2]
            ) / (2 * c4 * rho * lp(U1(rho) + 2 * rho * U2(rho)))
        return nothing
    end
end

function dv1(; λ, v0, U1, rhomin, d=3, n=1, η=0)
    c1 = (2 + d - η) / (pi^(d / 2) * (d * (2 + d) * gamma(d / 2) * 2^(d - 1)))
    return v0 * (c1^-1) * (λ - d) * (1 + U1(rhomin))^2 / n
end

function dv1_CS(; λ, v0, U1, rhomin, d=3, n=1, η=0)
    c3 = (2 - η) * gamma(2 - d / 2) / (2 * (4 * π)^(d / 2))
    return (v0 * (-d + λ) * (abs(1 + U1(rhomin)))^(2 - d / 2)) / (c3 * n)
end

# function dv1_r(; λ, v0, U1, rhomin, d=3, n=1, η=0)
#     c4 = (1 / 2) * (1 / (4 * pi)^(d / 2)) * (1 / gamma(d / 2))
#     return ((d - λ) * v0) / (c4 * n * lp(U1(rhomin)))
# end

function dv1_rfun(d, n, η)
    lp = lpfun(d, n, η)
    return function dv1_r(; λ, v0, U1, rhomin, d=3, n=1, η=0)
        c4 = (1 / 2) * (1 / (4 * pi)^(d / 2)) * (1 / gamma(d / 2))
        return ((d - λ) * v0) / (c4 * n * lp(U1(rhomin)))
    end
end

function etafun(usol, d)
    rho0 = find_zero(x -> usol.sol(x)[2], (0.1, 5.0))
    upp = usol.sol(rho0)[3]
    vd = 1 / (2^(d + 1) * pi^(d / 2) * gamma(d / 2))
    return 16 * vd / d * (rho0 * upp^2) / (1 + 2 * rho0 * upp)^2
end

function etafun_CS(usol, d)
    rho0 = find_zero(x -> usol.sol(x)[2], (0.1, 5.0))
    upp = usol.sol(rho0)[3]
    c = rho0 * upp
    if d ≈ 2.0
        return 2 / (
            1 +
            (8 * c * (1 + 2 * c) * pi * rho0) /
            (2 * c * (1 + c) - (1 + 2 * c) * log(1 + 2 * c))
        )
    else
        return (
            1 / 2 -
            (
                2^(1 + d) *
                c *
                (1 + 2 * c)^2 *
                pi^(-1 + d / 2) *
                rho0 *
                gamma(d / 2) *
                sin((d * pi) / 2)
            ) / (
                (1 + 2 * c)^(d / 2) * (-2 + c * (-6 + d)) +
                (1 + 2 * c)^2 * (2 + c * (-2 + d))
            )
        )^(-1)
    end
end

function etafunsigma(usol, d)
    rho0 = find_zero(x -> usol.sol(x)[2], (0.00001, 5.0))
    upp = usol.sol(rho0)[3]
    uppp = usol.sol(rho0, Val{1})[3]
    vd = 1 / (2^(d + 1) * pi^(d / 2) * gamma(d / 2))
    return 16 * vd / d * (3 * rho0 * upp^2 + 2 * rho0^2 * uppp)^2 / (1 + 2 * rho0 * upp)^4
end

function etafun_r(usol, r, rp, rpp)
    rho0 = find_zero(x -> usol.sol(x)[2], (0.1, 5.0))
    upp = usol.sol(rho0)[3]
    d = usol.d
    B1B2 = quadgk(x -> B1B2fun(d, rho0, upp, r, rp, rpp, x), 0.01, 0.1, 1.0, 10.0; order=10)[1]
    A1A2 = quadgk(x -> A1A2fun(d, rho0, upp, r, rp, rpp, x), 0.01, 0.1, 1.0, 10.0; order=10)[1]
    return -B1B2/(1+A1A2)
end
