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

x̂(i::Integer) = 3 * i / n

function e(i::Integer, x::Real)::Real
    if (x̂(i - 1) < x < x̂(i))
        return (x - x̂(i - 1)) / (3 / n)
    elseif (x̂(i) < x < x̂(i + 1))
        return (x̂(i + 1) - x) / (3 / n)
    end
    return 0
end

function e′(i::Integer, x::Real)::Real
    if (x̂(i - 1) < x < x̂(i))
        return 1 / (3 / n)
    elseif (x̂(i) < x < x̂(i + 1))
        return -1 / (3 / n)
    end
    return 0
end

" returns B(e_i, e_j) "
function B(i::Integer, j::Integer)::Real
    i, j = minmax(i, j)
    if (j == i)
        return -(∫(x -> e′(i, x)^2, x̂(i - 1), x̂(i)) +
                 ∫(x -> e′(i, x)^2, x̂(i), x̂(i + 1)))
    elseif (j == i + 1)
        return -(∫(x -> e′(i, x) * e′(j, x), x̂(i - 1), x̂(i)) +
                 ∫(x -> e′(i, x) * e′(j, x), x̂(i), x̂(j)) +
                 ∫(x -> e′(i, x) * e′(j, x), x̂(j), x̂(j + 1)))
    else
        return 0
    end
end

" returns L(e_i) "
function L(i::Integer)::Real
    result = ∫(x -> Φ̃′(x) * e′(i, x), x̂(i - 1), x̂(i)) +
             ∫(x -> Φ̃′(x) * e′(i, x), x̂(i), x̂(i + 1))
    if (x̂(i + 1) >= 1 && x̂(i - 1) <= 2)
        result += 4π * G * (∫(x -> e(i, x), x̂(i - 1), x̂(i)) +
                            ∫(x -> e(i, x), x̂(i), x̂(i + 1)))
    end
    return result
end

A = [B(i, j) for i in 1:n-1, j in 1:n-1]
Y = [L(i) for i in 1:n-1]
weights = A \ Y

function Φ(x::Real)
    return dot(weights, [e(i, x) for i in 1:n-1]) + Φ̃(x)
end

x = range(0, 3, 1000)
plot(x, Φ.(x), label=L"\Phi", lc=:black, lw=2, legend=false, xguidefontsize=14, yguidefontsize=14)
xlabel!(L"x")
ylabel!(L"\Phi(x)")
savefig("plot.pdf")
