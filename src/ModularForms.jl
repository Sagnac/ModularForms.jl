module ModularForms

using ComplexColor
import ComplexColor: GLMakie.Screen, complex_plot
using SpecialFunctions: zeta

export ModularForm, S, E, j, q, complex_plot

const ∑ = sum
const ∏ = prod
const ζ = zeta
const t = 275

ei2pi(x) = cispi(2x)

abstract type ModularForm end

struct E <: ModularForm
    k::Int
    c::Float64
    t::Int
    function E(k, t)
        (isodd(k) || k < 2) && error("Weight must be even and greater than zero.")
        c = 2 / ζ(1 - k)
        new(k, c, t)
    end
end

E(k) = E(k, t)

struct ModularFunction{F} <: ModularForm
    f::F
end

(F::ModularFunction)(q) = F.f(q)

S(a, q, t) = ∑(n ^ a * q ^ n / (1 - q ^ n) for n ∈ 1:t)

function (E_k::E)(q)
    (; k, c, t) = E_k
    1 + c * S(k - 1, q, t)
end

E4 = E(4)

Δ(q, t = t) = q * ∏((1 - q ^ n) ^ 24 for n ∈ 1:t)

g2(q) = (4 * π^4 / 3) * E4(q)

g3(q) = (8 * π^6 / 27) * E(6)(q)

j_(q, t = t) = E(4, t)(q) ^ 3 / Δ(q, t)

# function j_(q, t = t)
    # g2_3 = g2(q, t) ^ 3
    # 1728 * g2_3 / (g2_3 - 27 * g3(q, t) ^ 2)
# end

j = ModularFunction(j_)

q(τ) = ei2pi(τ)

function logclamp!(w)
    max_w = maximum(filter(isfinite, w))
    map!(w, w) do x
        log1p(isnan(x) ? max_w : clamp(x, 0, max_w))
    end
end

function complex_plot(f::ModularForm, n = 500)

    x = y = range(-1, 1, n)
    q = complex.(x, y')
    M = f.(q)
    M[abs.(q) .≥ 1] .= 0
    u, v = reim(M)
    logclamp!(u)
    logclamp!(v)
    figure = ComplexColor.Figure(; size = (2n, n))
    titlesize = 21
    u_axis = ComplexColor.Axis(figure[1,1]; title = L"Re_+(f(q))", titlesize)
    v_axis = ComplexColor.Axis(figure[1,2]; title = L"Im_+(f(q))", titlesize)
    ComplexColor.hidedecorations!(u_axis)
    ComplexColor.hidedecorations!(v_axis)
    kwargs = (; interpolate = false, colormap = :pink)

    ComplexColor.image!(u_axis, u; kwargs...)
    ComplexColor.image!(v_axis, v; kwargs...)

    display(Screen(), figure)

end

end
