module ModularForms

using ComplexColor
import ComplexColor: GLMakie.Screen, complex_plot
using SpecialFunctions: zeta

export ModularForm, ModularFunction, S, G, E, j, q, complex_plot

const ∑ = sum
const ∏ = prod
const ζ = zeta
const t = 275

ei2pi(x) = cispi(2x)

function γ(k)
    (isodd(k) || k < 4) && throw(DomainError(k, "\nrequired: k ≥ 4, k even\n"))
    2 / ζ(1 - k)
end

abstract type ModularForm end

# Eisenstein series
struct G <: ModularForm
    k::Int
    b::Float64
    c::Float64
    t::Int
    function G(k::Int, t::Int)
        b = 2ζ(k)
        c = b * γ(k) # == 2(2πi)^k / (k - 1)!
                     # == 2.0^(k+1) * π^k * (-1.0)^(k/2) / ∏(2.0:float(k-1))
        new(k, b, c, t)
    end
end

# normalized Eisenstein series
struct E <: ModularForm
    k::Int
    c::Float64
    t::Int
    function E(k::Int, t::Int)
        c = γ(k)
        new(k, c, t)
    end
end

G(k::Int) = G(k, t)

E(k::Int) = E(k, t)

struct ModularFunction{F} <: ModularForm
    f::F
end

(F::ModularFunction)(q) = F.f(q)

function S(a, q, t)
    s = 0.0im
    for n ∈ 1:t
        qⁿ = q^n
        s += n^a * qⁿ / (1 - qⁿ)
    end
    return s
end

# Fourier series expansion of the Eisenstein series
function eisenstein(b, c, f, q)
    (; k, t) = f
    b + c * S(k - 1, q, t)
end

(Gₖ::G)(q) = eisenstein(Gₖ.b, Gₖ.c, Gₖ, q)

(Eₖ::E)(q) = eisenstein(1, Eₖ.c, Eₖ, q)

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

    figure

end

end
