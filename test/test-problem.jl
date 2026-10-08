@testitem "PadeProblem degrees" tags = [:unit] begin
    c = collect(1.0:7.0)

    prob = PadeProblem(c)
    @test (prob.m, prob.n) == (3, 3)
    @test PadeProblem(c[1:6]).m == 3
    @test PadeProblem(c[1:6]).n == 2

    @test (PadeProblem(c; m = 4).m, PadeProblem(c; m = 4).n) == (4, 2)
    @test (PadeProblem(c; n = 1).m, PadeProblem(c; n = 1).n) == (5, 1)
    @test (PadeProblem(c; m = 2, n = 2).m, PadeProblem(c; m = 2, n = 2).n) == (2, 2)
    @test PadeProblem(c, 2, 2).c === c

    @test PadeProblem([1 // 1, 1 // 2, 1 // 6]) isa PadeProblem{Rational{Int}}
    @test PadeProblem(view(c, 1:5)).m == 2
    @test occursin("type (3, 3)", sprint(show, prob))
end

@testitem "PadeProblem validation" tags = [:unit] begin
    c = collect(1.0:5.0)

    @test_throws ArgumentError PadeProblem(c; m = 3, n = 2)
    @test_throws ArgumentError PadeProblem(c; m = 6)
    @test_throws ArgumentError PadeProblem(c; m = -1, n = 1)
    @test_throws ArgumentError PadeProblem(Float64[])
    @test PadeProblem(c; m = 2, n = 2) isa PadeProblem
end

@testitem "Return codes" tags = [:unit] begin
    @test successful_retcode(ReturnCode.Success)
    @test successful_retcode(ReturnCode.DegreeReduced)
    @test !successful_retcode(ReturnCode.Breakdown)
end
