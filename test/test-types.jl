@testitem "RationalApproximant" tags = [:unit] begin
    # (1 + z/2) / (1 - z/2), the type (1, 1) approximant of exp
    r = RationalApproximant([1.0, 0.5], [1.0, -0.5], 1, 1, ReturnCode.Success)

    @test degrees(r) == (1, 1)
    @test r(0.0) == 1.0
    @test r(1.0) == 3.0
    @test r(1.0im) ≈ (1 + 0.5im) / (1 - 0.5im)
    @test successful_retcode(r)
    @test occursin("type (1, 1)", sprint(show, r))

    reduced = RationalApproximant([2.0], [1.0], 3, 3, ReturnCode.DegreeReduced)
    @test degrees(reduced) == (0, 0)
    @test (reduced.m, reduced.n) == (3, 3)
    @test reduced(10.0) == 2.0

    exact = RationalApproximant([1 // 1, 1 // 2], [1 // 1, -1 // 2], 1, 1, ReturnCode.Success)
    @test exact(1 // 2) === 5 // 3

    @test_throws ArgumentError RationalApproximant([1.0], [2.0], 0, 0, ReturnCode.Success)
    @test_throws ArgumentError RationalApproximant([1.0], Float64[], 0, 0, ReturnCode.Success)
end

@testitem "ContinuedFraction" tags = [:unit] begin
    # 1 / (1 - z / (1 + z/2)) = (1 + z/2) / (1 - z/2)
    cf = ContinuedFraction(Float64[], [1.0, -1.0, 0.5], 1, 1, ReturnCode.Success)

    @test degrees(cf) == (1, 1)
    @test cf(0.0) == 1.0
    @test cf(1.0) ≈ 3.0
    @test cf(0.3 + 0.2im) ≈ (1 + (0.3 + 0.2im) / 2) / (1 - (0.3 + 0.2im) / 2)
    @test occursin("type (1, 1)", sprint(show, cf))

    # Number of partial numerators against degrees: staircase [0/0], [0/1], [1/1], [1/2]
    for (L, type) in ((1, (0, 0)), (2, (0, 1)), (3, (1, 1)), (4, (1, 2)))
        @test degrees(ContinuedFraction(Float64[], ones(L), 0, 0, ReturnCode.Success)) == type
    end

    # 2 + 3z + z² / (1 + z): a prefix shifts the fraction by its length
    shifted = ContinuedFraction([2.0, 3.0], [1.0, 1.0], 2, 1, ReturnCode.Success)
    @test degrees(shifted) == (2, 1)
    @test shifted(0.5) ≈ 2 + 3 * 0.5 + 0.5^2 / 1.5

    polynomial = ContinuedFraction([1.0, 2.0], Float64[], 1, 0, ReturnCode.Success)
    @test degrees(polynomial) == (1, 0)
    @test polynomial(2.0) == 5.0

    exact = ContinuedFraction(Rational{Int}[], [1 // 1, -1 // 1, 1 // 2], 1, 1, ReturnCode.Success)
    @test exact(1 // 2) === 5 // 3
end

@testitem "PointEstimate" tags = [:unit] begin
    est = PointEstimate(2.718, 1.0, 1.0e-3, 2, 2, ReturnCode.Success)

    @test value(est) == 2.718
    @test degrees(est) == (2, 2)
    @test successful_retcode(est)
    @test occursin("2.718", sprint(show, est))
end
