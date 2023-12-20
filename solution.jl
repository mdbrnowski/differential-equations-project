using SparseArrays
using LinearAlgebra
using LaTeXStrings
using Plots
include("integration.jl")

const n::Integer = 100
# const G::Real = 6.6743015e-11
const G::Real = 1

Φ̃(x::Real) = 5 - x / 3

Φ̃′(x::Real) = -1 / 3

xi(i::Integer) = 3 * i / n

function e(i::Integer, x::Real)::Real
    if (xi(i - 1) < x < xi(i))
        return (x - xi(i - 1)) / (3 / n)
    elseif (xi(i) < x < xi(i + 1))
        return (xi(i + 1) - x) / (3 / n)
    end
    return 0
end

function e′(i::Integer, x::Real)::Real
    if (xi(i - 1) < x < xi(i))
        return 1 / (3 / n)
    elseif (xi(i) < x < xi(i + 1))
        return -1 / (3 / n)
    end
    return 0
end

" returns B(e_i, e_j) "
function B(i::Integer, j::Integer)
    i, j = minmax(i, j)
    if (j == i)
        return -(∫(x -> e′(i, x)^2, xi(i - 1), xi(i)) +
                 ∫(x -> e′(i, x)^2, xi(i), xi(i + 1)))
    elseif (j == i + 1)
        return -(∫(x -> e′(i, x) * e′(j, x), xi(i - 1), xi(i)) +
                 ∫(x -> e′(i, x) * e′(j, x), xi(i), xi(j)) +
                 ∫(x -> e′(i, x) * e′(j, x), xi(j), xi(j + 1)))
    else
        return 0
    end
end

" returns L(e_i) "
function L(i::Integer)::Real
    result = ∫(x -> Φ̃′(x) * e′(i, x), xi(i - 1), xi(i)) +
             ∫(x -> Φ̃′(x) * e′(i, x), xi(i), xi(i + 1))
    if (xi(i + 1) >= 1 && xi(i - 1) <= 2)
        result += 4π * G * (∫(x -> e(i, x), xi(i - 1), xi(i)) +
                            ∫(x -> e(i, x), xi(i), xi(i + 1)))
    end
    return result
end

A = zeros(Real, (n - 1, n - 1))
for i in 1:n-1, j in 1:n-1
    A[i, j] = B(i, j)
end

Y = zeros(Real, (n - 1, 1))
for i in 1:n-1
    Y[i, 1] = L(i)
end

weights = A\Y

function Φ(x::Real)
    return dot(weights, [e(i, x) for i in 1:n-1]) + Φ̃(x)
end

x = range(0, 3, 100)
plot(x, Φ.(x), label=L"\Phi", lc=:black, lw=2, legend=false, xguidefontsize=14, yguidefontsize=14)
xlabel!(L"x")
ylabel!(L"\Phi(x)")
savefig("plot.pdf")
