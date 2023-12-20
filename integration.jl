using LinearAlgebra

function gauss_legendre_quadrature(f::Function)::Real
    xs = [-1 / sqrt(3), 1 / sqrt(3)]
    weights = [1, 1]
    return dot(weights, f.(xs))
end

" Integrates f on range (a, b). Uses Gauss-Legendre quadrature with n = 3. "
function integrate(f::Function, a::Real, b::Real)::Real
    g(x) = f((x + 1) / 2 * (b - a) + a)
    result = gauss_legendre_quadrature(g)
    return result * (b - a) / 2
end

const ∫ = integrate