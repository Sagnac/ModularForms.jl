module ModularForms

using ComplexColor
using ComplexColor: GLMakie.Screen
using Main: ei2pi

export S, g2, g3, j, q, figure

const ∑ = sum
const ∏ = prod

# function d(n)
    # v = [1]
    # f = 2
    # while f < n
        # iszero(n % f) && push!(v, f)
        # f += 1
    # end
    # push!(v, n)
    # return v
# end

# σ(a, n) = sum(d(n) .^ a)

t = 275

function S(a, q, t = t)
    ∑(n ^ a * q ^ n / (1 - q ^ n) for n = 1:t)
    # sum(σ(a, n) * q ^ n for n = 1:t)
end

E4(q, t = t) = 1 + 240 * S(3, q, t)

Δ(q, t = t) = q * ∏((1 - q ^ n) ^ 24 for n = 1:t)

g2(q, t = t) = (4 * π^4 / 3) * E4(q, t)

g3(q, t = t) = (8 * π^6 / 27) * (1 - 504 * S(5, q, t))

j(q, t = t) = E4(q, t) ^ 3 / Δ(q, t)

# function j(q, t = t)
    # g2_3 = g2(q, t) ^ 3
    # 1728 * g2_3 / (g2_3 - 27 * g3(q, t) ^ 2)
# end

q(τ) = ei2pi(τ)

len = 500

x = y = range(-1, 1, len)

_q = complex.(x, y')

mask = abs.(_q) .≥ 1

_g2 = g2.(_q)
# _g3 = g3.(_q)
# _j = j.(_q)

_g2[mask] .= 0
# _g3[mask] .= 0
# _j[mask] .= 0

Re_g2 = real(_g2)
# Re_g3 = real(_g3)
# Re_j = real(_j)

Im_g2 = imag(_g2)

clamp!(Re_g2, 0, maximum(Re_g2))
# clamp!(Re_g3, 0, maximum(Re_g3))
# clamp!(Re_j, 0, maximum(isfinite, Re_j))

clamp!(Im_g2, 0, maximum(Im_g2))

map!(x -> log(x + 1), Re_g2, Re_g2)
# map!(x -> log(x + 1), Re_g3, Re_g3)
# map!(x -> log(x + 1), Re_j, Re_j)

map!(x -> log(x + 1), Im_g2, Im_g2)

figure = ComplexColor.Figure(; size = (2len, len))
titlesize = 21
axis_Re = ComplexColor.Axis(figure[1,1]; title = L"Re(g_2(q))", titlesize)
axis_Im = ComplexColor.Axis(figure[1,2]; title = L"Im(g_2(q))", titlesize)
ComplexColor.hidedecorations!(axis_Re)
ComplexColor.hidedecorations!(axis_Im)
kwargs = (; interpolate = false, colormap = :pink)

ComplexColor.image!(axis_Re, Re_g2; kwargs...)
ComplexColor.image!(axis_Im, Im_g2; kwargs...)
# ComplexColor.image(Im_g2)
# ComplexColor.image(angle.(_g2))

# display(figure)

end
