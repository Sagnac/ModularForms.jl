using ComplexColor
using ComplexColor: GLMakie.Screen

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
    sum(n ^ a * q ^ n / (1 - q ^ n) for n = 1:t)
    # sum(σ(a, n) * q ^ n for n = 1:t)
end

g2(q, t = t) = (4 * π^4 / 3) * (1 + 240 * S(3, q, t))

g3(q, t = t) = (8 * π^6 / 27) * (1 - 504 * S(5, q, t))

function j(q, t = t)
    g2_3 = g2(q, t) ^ 3
    1728 * g2_3 / (g2_3 - 27 * g3(q, t) ^ 2)
end

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

image_1 = ComplexColor.image(Re_g2, interpolate = false)
# ComplexColor.image(Im_g2, interpolate = false)
# ComplexColor.image(angle.(_g2), interpolate = false)
image_2 = ComplexColor.image(Im_g2, interpolate = false)

display(Screen(), image_1)
display(Screen(), image_2)
