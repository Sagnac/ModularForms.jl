module ModularForms

using ComplexColor
import ComplexColor: GLMakie.Screen, complex_plot, ComplexArray
using SpecialFunctions: zeta

export ModularForm, ModularFunction, S, G, E, j, q, ReIm, ModArg

for name in names(ComplexColor)
    @eval export $name
end

const ∑ = sum
const ∏ = prod
const ζ = zeta
const t = Ref(275)
const n = Ref(500)

q(τ) = cispi(2τ)

function γ(k)
    (isodd(k) || k < 4) && throw(DomainError(k, "\nrequired: k ≥ 4, k even\n"))
    2 / ζ(1 - k)
end

abstract type ModularForm end

abstract type EisensteinSeries <: ModularForm end

struct ReIm end
struct ModArg end

# Eisenstein series
struct G <: EisensteinSeries
    k::Int
    b::Float64
    c::Float64
    function G(k::Int)
        b = 2ζ(k)
        c = b * γ(k) # == 2(2πi)^k / (k - 1)!
                     # == 2.0^(k+1) * π^k * (-1.0)^(k/2) / ∏(2.0:float(k-1))
        new(k, b, c)
    end
end

# normalized Eisenstein series
struct E <: EisensteinSeries
    k::Int
    c::Float64
    function E(k::Int)
        c = γ(k)
        new(k, c)
    end
end

struct ModularFunction{F} <: ModularForm
    f::F
end

(F::ModularFunction)(q::Number) = F.f(q)
(F::ModularFunction)(q::Number, t::Int) = F.f(q, t)

# Lambert series
function S(a::Number, q::Number, t::Int)
    s = 0.0im
    for n ∈ 1:t
        qⁿ = q^n
        s += n^a * qⁿ / (1 - qⁿ)
    end
    return s
end

# Fourier series expansion of the Eisenstein series
function (f::EisensteinSeries)(q::Number, t::Int = t[])
    (; b, c, k) = f
    b + c * S(k - 1, q, t)
end

Base.getproperty(f::E, name::Symbol) = name == :b ? 1 : getfield(f, name)

E4 = E(4)
E6 = E(6)

# normalized modular discriminant η²⁴(q) [regular discriminant divided by (2π)¹²]
# where η is the Dedekind eta function
Δ(q, t = t[]) = q * ∏((1 - q ^ n) ^ 24 for n ∈ 1:t)

g2(q, t = t[]) = (4 * π^4 / 3) * E4(q, t)

g3(q, t = t[]) = (8 * π^6 / 27) * E6(q, t)

# j(q) = 12^3 * J(q) where J is Felix Klein's Absolute Invariant
j_invariant(q::Number, t::Int = t[]) = E4(q, t) ^ 3 / Δ(q, t)

# function j_invariant(q, t = t[])
    # g2_3 = g2(q, t) ^ 3
    # 1728 * g2_3 / (g2_3 - 27 * g3(q, t) ^ 2)
# end

j = ModularFunction(j_invariant)

set_title(f::T) where T <: EisensteinSeries = string(T, "_", f.k)
set_title(f::ModularForm) = f ≡ j ? "j" : "f"

function logclamp!(w)
    max_w = maximum(filter(isfinite, w))
    map!(w, w) do x
        log1p(isnan(x) ? max_w : clamp(x, 0, max_w))
    end
end

function complex_plot(f::ModularForm, ::Type{ReIm}, t::Int = t[], n::Int = n[])
    x = y = range(-1, 1, n)
    q = complex.(x, y')
    M = f.(q, t)
    M[abs.(q) .≥ 1] .= 0
    u, v = reim(M)
    logclamp!(u)
    logclamp!(v)
    figure = ComplexColor.Figure(; size = (2n, n))
    titlesize = 21
    f_title = set_title(f)
    u_axis = ComplexColor.Axis(figure[1,1]; title = L"Re_+(%$f_title(q))", titlesize)
    v_axis = ComplexColor.Axis(figure[1,2]; title = L"Im_+(%$f_title(q))", titlesize)
    ComplexColor.hidedecorations!(u_axis)
    ComplexColor.hidedecorations!(v_axis)
    kwargs = (; interpolate = false, colormap = :pink)
    ComplexColor.image!(u_axis, u; kwargs...)
    ComplexColor.image!(v_axis, v; kwargs...)
    figure
end

function complex_plot(f::ModularForm, t::Int = t[], n::Int = n[])
    complex_plot(f, ReIm, t, n)
end

function complex_plot(f::ModularForm, ::Type{ModArg},
                      color::ComplexColor.Spaces = ComplexColor.default,
                      t::Int = t[], n::Int = n[]; title = L"%$(set_title(f))(q)")
    x = y = range(-1, 1, n)
    complex_plot(x, y, q -> abs(q) < 1 ? f(q, t) : 0.0im, color; title)
end

# f(τ) methods

function complex_plot(x::AbstractVector, y::AbstractVector, f::ModularForm,
                      color::ComplexColor.Spaces = ComplexColor.default;
                      title = L"%$(set_title(f))(τ)")
    complex_plot(x, y, f ∘ q, color; title)
end

function complex_plot(x::Interval, y::Interval, f::ModularForm,
                      color::ComplexColor.Spaces = ComplexColor.default;
                      title = L"%$(set_title(f))(τ)")
    complex_plot(x, y, f ∘ q, color; title)
end

function complex_plot(x::AbstractVector, y::Interval, f::ModularForm,
                      color::ComplexColor.Spaces = ComplexColor.default;
                      title = L"%$(set_title(f))(τ)")
    complex_plot(x, y, f ∘ q, color; title)
end

function complex_plot(x::Interval, y::AbstractVector, f::ModularForm,
                      color::ComplexColor.Spaces = ComplexColor.default;
                      title = L"%$(set_title(f))(τ)")
    complex_plot(x, y, f ∘ q, color; title)
end

function complex_plot(z::ComplexArray, f::ModularForm,
                      color::ComplexColor.Spaces = ComplexColor.default;
                      title = L"%$(set_title(f))(τ)")
    complex_plot(z, f ∘ q, color; title)
end

end
