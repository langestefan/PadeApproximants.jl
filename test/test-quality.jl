@testitem "Aqua quality assurance" tags = [:quality] begin
    using Aqua: Aqua
    Aqua.test_all(PadeApproximants)
end

@testitem "JET static analysis" tags = [:quality] begin
    using JET: JET
    JET.test_package(PadeApproximants; target_modules = (PadeApproximants,))
end
