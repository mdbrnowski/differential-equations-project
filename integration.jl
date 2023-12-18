using LinearAlgebra

function gauss_legendre_quadrature(f::Function, n::Int64)::Float64
    if n == 2
        xs = [-1 / sqrt(3), 1 / sqrt(3)]
        weights = [1, 1]
    elseif n == 3
        xs = [-sqrt(3 / 5), 0, sqrt(3 / 5)]
        weights = [5 / 9, 8 / 9, 5 / 9]
    else
        error("n should be equal to 2 or 3")
    end
    result = dot(weights, f.(xs))
    return result
end

" Integrates f on range (a, b). Uses Gauss-Legendre quadrature with n = 3. "
function integrate(f::Function, a::Float64, b::Float64)::Float64
    g(x) = f((x + 1) / 2 * (b - a) + a)
    result = gauss_legendre_quadrature(g, 3)
    return result * (b - a) / 2
end