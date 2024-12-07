using ChemicalKinetics
using Test

function test_greet()
    @test 1 == 1
end

@testset "ChemicalKinetics.jl" begin
	test_greet()
end