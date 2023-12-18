using LinearAlgebra

function gauss_legendre_quadrature_on_standard_interval(f, n)
    if n == 2
        println(2)
        xs = [-1/sqrt(3), 1/sqrt(3)]
        weights = [1, 1]
    elseif n == 3
        println(3)
        xs = [-sqrt(3/5), 0, sqrt(3/5)]
        weights = [5/9, 8/9, 5/9]
    else
        error("n should be equal to 2 or 3")
    end
    result = dot(weights, map(f, xs))
    return result
end

function gauss_legendre_quadrature(f, a, b)
    g(x) = f((x + 1) / 2 * (b - a) + a)
    result = gauss_legendre_quadrature_on_standard_interval(g, 3)
    return result * (b - a) / 2
end