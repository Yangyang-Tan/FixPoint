using QuadGK, Interpolations, Dierckx
function r_exp(x, n)
    if x^n > 300.0
        return 0 * x
    elseif 0 <= x^n < 0.0001
        return 1 / x - x^(n - 1) / 2
    else
        x^(-1 + n) / (-1 + exp(x^n))
    end
end

function rp_exp(x, n)
    if x^n > 300.0
        return 0 * x
    elseif 0 <= x^n < 0.0001
        return -x^(-2) - ((-1 + n) * x^(-2 + n)) / 2 + ((-1 + 2 * n) * x^(-2 + 2 * n)) / 12
    else
        ((-1 + n) * x^(-2 + n)) / (-1 + exp(x^n)) -
        (n * x^(-2 + 2 * n) * exp(x^n)) / (-1 + exp(x^n))^2
    end
end

function rpp_exp(x, n)
    if x^n > 300.0
        return 0 * x
    elseif 0 <= x^n < 0.0001
        return 2 / x^3 - ((-2 + n) * (-1 + n) * x^(-3 + n)) / 2 +
               ((-2 + 2 * n) * (-1 + 2 * n) * x^(-3 + 2 * n)) / 12 -
               ((-2 + 4 * n) * (-1 + 4 * n) * x^(-3 + 4 * n)) / 720
    else
        ((-2 + n) * (-1 + n) * x^(-3 + n)) / (-1 + exp(x^n)) -
        (2 * (-1 + n) * n * x^(-3 + 2 * n) * exp(x^n)) / (-1 + exp(x^n))^2 +
        x^(-1 + n) * (
            -(((-1 + n) * n * x^(-2 + n) * exp(x^n)) / (-1 + exp(x^n))^2) -
            (n^2 * x^(-2 + 2 * n) * exp(x^n)) / (-1 + exp(x^n))^2 +
            (2 * n^2 * x^(-2 + 2 * n) * exp(2 * x^n)) / (-1 + exp(x^n))^3
        )
    end
end

function r_exp(n)
    return x -> r_exp(x, n)
end

function rp_exp(n)
    return x -> rp_exp(x, n)
end

function rpp_exp(n)
    return x -> rpp_exp(x, n)
end

function lp_ker(y, d, n, η, ωp1)
    r = r_exp
    rp = rp_exp
    return (y^(d / 2) * (η * r(y, n) + 2 * y * rp(y, n))) / (-1 + y + ωp1 + y * r(y, n))^2
end

function lpfun_val(d, n, η, ωp1)
    return quadgk(
        y -> lp_ker(y, d, n, η, ωp1),
        0.0,
        Inf;
        # maxevals=10^8,
        order=10,
        # atol=1e-10,
        # rtol=1e-10,
    )[1]
end

function lpfun(d, n, η)
    T = eltype([d, n, η])
    lprang = exp.(range(T(log(10^-6)); stop=T(log(10^12)), length=2000))
    lprang_itp = range(log(10^-6); stop=log(10^12), length=2000)
    itp_scale = scale(
        interpolate(lpfun_val.(d, n, η, lprang), BSpline(Linear())), lprang_itp
    )
    return x -> itp_scale(log(x + 1))
end



function lpfun(d, n, η)
    T = eltype([d, n, η])
    lprang = exp.(range(T(log(10^-4)); stop=T(log(10^4)), length=1000))
    lprang_itp = range(log(10^-4); stop=log(10^4), length=1000)
    itp_scale = scale(
        interpolate(lpfun_val.(d, n, η, lprang), BSpline(Linear())), lprang_itp
    )
    lprang2 = exp.(range(T(log(10^-8)); stop=T(log(10^-4)), length=50))
    lpval2 = -lpfun_val.(d, n, η, lprang2)
    lprang3 = exp.(range(T(log(10^4)); stop=T(log(10^8)), length=50))
    lpval3 = -lpfun_val.(d, n, η, lprang3)
    fun2=curve_fit(PowerFit,lprang2, lpval2)
    fun3=curve_fit(PowerFit,lprang3, lpval3)
    return x -> begin
        if x+1 < 10^-4
            -fun2(x+1)
        elseif x+1 > 10^4
            -fun3(x+1)
        else
            itp_scale(log(x + 1))
        end
    end
end


# function lpfun(d, n, η)
#     T = eltype([d, n, η])
#     lprang = exp.(range(T(log(10^-6)); stop=log(T(log(10^10))), length=1000))
#     lprang_itp = range(T(log(10^-6)); stop=log(T(log(10^10))), length=1000)
#     itp_scale = Spline1D(lprang_itp,lpfun_val.(d, n, η, lprang))
#     return x -> itp_scale(log(x + 1))
# end

function B1B2fun(d, ρbar, upp, r::Function, rp::Function, rpp::Function, x)
    return (
        2^(2 - d) *
        upp^2 *
        x^(-2 + d / 2) *
        ρbar *
        rp(x) *
        (
            x *
            (x + 2 * upp * ρbar + x * r(x)) *
            (
                d * (1 + r(x)) * (1 + r(x) + x * rp(x)) - 4 * (1 + r(x) + x * rp(x))^2 +
                2 * x * (1 + r(x)) * (2 * rp(x) + x * rpp(x))
            ) +
            x *
            (1 + r(x)) *
            (
                d * (x + 2 * upp * ρbar + x * r(x)) * (1 + r(x) + x * rp(x)) -
                4 * x * (1 + r(x) + x * rp(x))^2 +
                2 * x * (x + 2 * upp * ρbar + x * r(x)) * (2 * rp(x) + x * rpp(x))
            )
        )
    ) / (
        d * pi^(d / 2) * gamma(d / 2) * (1 + r(x))^3 * (x + 2 * upp * ρbar + x * r(x))^3
    )
end

function A1A2fun(d, ρbar, upp, r::Function, rp::Function, rpp::Function, x)
    return (
        2^(1 - d) *
        upp^2 *
        x^(-3 + d / 2) *
        ρbar *
        r(x) *
        (
            x *
            (x + 2 * upp * ρbar + x * r(x)) *
            (
                d * (1 + r(x)) * (1 + r(x) + x * rp(x)) - 4 * (1 + r(x) + x * rp(x))^2 +
                2 * x * (1 + r(x)) * (2 * rp(x) + x * rpp(x))
            ) +
            x *
            (1 + r(x)) *
            (
                d * (x + 2 * upp * ρbar + x * r(x)) * (1 + r(x) + x * rp(x)) -
                4 * x * (1 + r(x) + x * rp(x))^2 +
                2 * x * (x + 2 * upp * ρbar + x * r(x)) * (2 * rp(x) + x * rpp(x))
            )
        )
    ) / (
        d * pi^(d / 2) * gamma(d / 2) * (1 + r(x))^3 * (x + 2 * upp * ρbar + x * r(x))^3
    )
end
