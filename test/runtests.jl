using ChemicalKinetics
using Test

@testset "Reader.jl" begin
    species = ChemicalKinetics.read_species("../data/CO2.json")

    @test "CO2" == species.name
    @test 7.31e-26 == species.mass
    @test 2 == species.dof_rot
    @test 4 == species.dof_vib
    @test 5 == species.Zrot

    @test 5.62e-10 == species.vhs.dref
    @test 273.0 == species.vhs.Tref
    @test 0.8 == species.vhs.omega

    @test 1 == species.vhs.alpha
    @test 18.1 == species.vhs.Zrotinf
    @test 91.5 == species.vhs.Tstar
    @test 25.48 == species.vhs.C1
    @test 177.98 == species.vhs.C2
end