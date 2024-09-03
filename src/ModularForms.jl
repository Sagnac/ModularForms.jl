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

function complex_plot(f::ModularForm, n = 500)

    x = y = range(-1, 1, n)
    q = complex.(x, y')
    M = f.(q)
    M[abs.(q) .≥ 1] .= 0
    Re, Im = reim(M)
    max_Re = maximum(filter(isfinite, Re))
    max_Im = maximum(filter(isfinite, Im))
    map!(x -> !isfinite(x) ? max_Re : clamp(x, 0, max_Re), Re, Re)
    map!(x -> !isfinite(x) ? max_Im : clamp(x, 0, max_Im), Im, Im)
    map!(x -> log(x + 1), Re, Re)
    map!(x -> log(x + 1), Im, Im)

    figure = ComplexColor.Figure(; size = (2n, n))
    titlesize = 21
    axis_Re = ComplexColor.Axis(figure[1,1]; title = L"Re_+(f(q))", titlesize)
    axis_Im = ComplexColor.Axis(figure[1,2]; title = L"Im_+(f(q))", titlesize)
    ComplexColor.hidedecorations!(axis_Re)
    ComplexColor.hidedecorations!(axis_Im)
    kwargs = (; interpolate = false, colormap = :pink)

    ComplexColor.image!(axis_Re, Re; kwargs...)
    ComplexColor.image!(axis_Im, Im; kwargs...)

    display(Screen(), figure)

end

end
