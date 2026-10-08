"""
    AbstractPadeResult

Supertype of the results returned by `solve`: [`RationalApproximant`](@ref),
[`ContinuedFraction`](@ref) and [`PointEstimate`](@ref). Every result has a `retcode` field
and stores the requested degrees in `m` and `n`.
"""
abstract type AbstractPadeResult end

successful_retcode(r::AbstractPadeResult) = successful_retcode(r.retcode)

"""
    degrees(r) -> (m, n)

Actual numerator and denominator degrees of the result `r`. These can be lower than the
requested degrees `(r.m, r.n)`.
"""
function degrees end

"""
    value(est::PointEstimate)

The estimated value of the series at `est.at`.
"""
function value end

"""
    RationalApproximant(p, q, m, n, retcode)

Rational function `p(z) / q(z)`, with coefficient vectors `p` and `q` stored constant term
first and `q[1] == 1`. `m` and `n` are the degrees that were requested; the actual degrees
are given by [`degrees`](@ref). Call the object to evaluate it: `r(z)`.
"""
struct RationalApproximant{T <: Number, V <: AbstractVector{T}} <: AbstractPadeResult
    p::V
    q::V
    m::Int
    n::Int
    retcode::ReturnCode.T

    function RationalApproximant(
            p::V,
            q::V,
            m::Integer,
            n::Integer,
            retcode::ReturnCode.T,
        ) where {T <: Number, V <: AbstractVector{T}}
        (isempty(p) || isempty(q)) &&
            throw(ArgumentError("numerator and denominator need at least one coefficient"))
        isone(first(q)) ||
            throw(ArgumentError("the denominator must have constant term 1"))
        return new{T, V}(p, q, m, n, retcode)
    end
end

degrees(r::RationalApproximant) = (length(r.p) - 1, length(r.q) - 1)

(r::RationalApproximant)(z::Number) = evalpoly(z, r.p) / evalpoly(z, r.q)

"""
    ContinuedFraction(prefix, a, m, n, retcode)

Approximant in continued-fraction form,

    prefix(z) + zᵏ · a₁ / (1 + a₂z / (1 + a₃z / (1 + …)))

where `prefix` holds the `k` coefficients of a polynomial of degree `k - 1`, constant term
first, and `a` the partial numerators. `m` and `n` are the degrees that were requested; the
actual degrees are given by [`degrees`](@ref). Call the object to evaluate it by backward
recurrence: `cf(z)`.
"""
struct ContinuedFraction{T <: Number, V <: AbstractVector{T}} <: AbstractPadeResult
    prefix::V
    a::V
    m::Int
    n::Int
    retcode::ReturnCode.T

    function ContinuedFraction(
            prefix::V,
            a::V,
            m::Integer,
            n::Integer,
            retcode::ReturnCode.T,
        ) where {T <: Number, V <: AbstractVector{T}}
        return new{T, V}(prefix, a, m, n, retcode)
    end
end

function degrees(cf::ContinuedFraction)
    k = length(cf.prefix)
    L = length(cf.a)
    L == 0 && return (max(k - 1, 0), 0)
    return (k + (L - 1) ÷ 2, L ÷ 2)
end

function (cf::ContinuedFraction{T})(z::Number) where {T}
    a = cf.a
    t = zero(z) / one(T)
    head = isempty(cf.prefix) ? t : evalpoly(z, cf.prefix) + t
    isempty(a) && return head
    for i in lastindex(a):-1:(firstindex(a) + 1)
        t = a[i] * z / (1 + t)
    end
    return head + z^length(cf.prefix) * (first(a) / (1 + t))
end

"""
    PointEstimate(value, at, abserr, m, n, retcode)

Value of the type `(m, n)` approximant at the single point `at`. `abserr` is the absolute
difference from the previous approximant in the same sequence, which serves as an error
estimate. Get the value with [`value`](@ref).
"""
struct PointEstimate{T <: Number, S <: Number, E <: Real} <: AbstractPadeResult
    value::T
    at::S
    abserr::E
    m::Int
    n::Int
    retcode::ReturnCode.T
end

degrees(est::PointEstimate) = (est.m, est.n)

value(est::PointEstimate) = est.value

function Base.show(io::IO, r::RationalApproximant{T}) where {T}
    print(io, "RationalApproximant{", T, "} of type ", degrees(r))
    return print(io, ", retcode: ", r.retcode)
end

function Base.show(io::IO, cf::ContinuedFraction{T}) where {T}
    print(io, "ContinuedFraction{", T, "} of type ", degrees(cf))
    return print(io, ", retcode: ", cf.retcode)
end

function Base.show(io::IO, est::PointEstimate)
    print(io, "PointEstimate(", est.value, " at ", est.at, ", abserr ≈ ", est.abserr)
    return print(io, ", retcode: ", est.retcode, ")")
end
