"""
    ReturnCode

Namespace for the return codes attached to every result, as `ReturnCode.T` values:

- `ReturnCode.Success`: the approximant of the requested type was computed.
- `ReturnCode.DegreeReduced`: the approximant is valid, but of lower type than requested,
  because the coefficients do not support the requested degrees.
- `ReturnCode.Breakdown`: the algorithm hit a singular system or a zero divisor. The result
  holds whatever was computed up to that point and should not be trusted.

See also [`successful_retcode`](@ref PadeApproximants.successful_retcode).
"""
module ReturnCode
    @enum T Success DegreeReduced Breakdown
end

"""
    successful_retcode(retcode::ReturnCode.T)
    successful_retcode(result)

Return `true` if the result can be used, that is, for `ReturnCode.Success` and
`ReturnCode.DegreeReduced`.
"""
function successful_retcode(retcode::ReturnCode.T)
    return retcode == ReturnCode.Success || retcode == ReturnCode.DegreeReduced
end

"""
    AbstractPadeAlgorithm

Supertype of the algorithms accepted by `solve(prob::PadeProblem, alg)`.
"""
abstract type AbstractPadeAlgorithm end

"""
    PadeProblem(c; m, n)
    PadeProblem(c, m, n)

Problem of finding the Padé approximant of type `(m, n)` to the power series with
coefficients `c`, stored constant term first (`c[1]` is c₀). `m` is the numerator degree
and `n` the denominator degree.

If neither degree is given, all coefficients are used for the (near-)diagonal approximant:
`n = (length(c) - 1) ÷ 2` and `m = length(c) - 1 - n`. If only one is given, the other takes
the remaining coefficients.

Throws an `ArgumentError` if `c` has fewer than `m + n + 1` coefficients. Coefficients
beyond the first `m + n + 1` are ignored.
"""
struct PadeProblem{T <: Number, V <: AbstractVector{T}}
    c::V
    m::Int
    n::Int

    function PadeProblem(c::AbstractVector{T}, m::Integer, n::Integer) where {T <: Number}
        Base.require_one_based_indexing(c)
        (m >= 0 && n >= 0) ||
            throw(ArgumentError("degrees must be non-negative, got m = $m, n = $n"))
        length(c) >= m + n + 1 || throw(
            ArgumentError(
                "a type ($m, $n) approximant needs $(m + n + 1) coefficients, " *
                    "got $(length(c))",
            ),
        )
        return new{T, typeof(c)}(c, m, n)
    end
end

function PadeProblem(
        c::AbstractVector{<:Number};
        m::Union{Nothing, Integer} = nothing,
        n::Union{Nothing, Integer} = nothing,
    )
    isempty(c) && throw(ArgumentError("the coefficient vector is empty"))
    N = length(c) - 1
    if m === nothing && n === nothing
        n = N ÷ 2
        m = N - n
    elseif m === nothing
        m = N - n
    elseif n === nothing
        n = N - m
    end
    return PadeProblem(c, m, n)
end

function Base.show(io::IO, prob::PadeProblem{T}) where {T}
    print(io, "PadeProblem{", T, "} of type (", prob.m, ", ", prob.n, ") with ")
    return print(io, length(prob.c), " coefficients")
end
