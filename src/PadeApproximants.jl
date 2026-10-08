module PadeApproximants

using CommonSolve: CommonSolve, init, solve, solve!
using LinearAlgebra: LinearAlgebra

export PadeProblem
export RationalApproximant, ContinuedFraction, PointEstimate
export ReturnCode, successful_retcode
export degrees, value
export init, solve, solve!

include("problem.jl")
include("types.jl")

end
