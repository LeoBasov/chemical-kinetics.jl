using ChemicalKinetics
using Test

@testset "ChemicalKinetics.jl" begin
    add_species!("../data/CO2.json")
    add_species!("../data/CH4.json")
    add_species!("../data/N2.json")

    @test 3 == length(ChemicalKinetics.state.species)
    
    @test 0.0 == ChemicalKinetics.state.mole_fractions["CO2"]
    @test 0.0 == ChemicalKinetics.state.mole_fractions["CH4"]
    @test 0.0 == ChemicalKinetics.state.mole_fractions["N2"]

    @test 3 == length(ChemicalKinetics.state.Tvib["CO2"])
    @test 4 == length(ChemicalKinetics.state.Tvib["CH4"])
    @test 1 == length(ChemicalKinetics.state.Tvib["N2"])
end

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

    @test 3 == length(species.vibmodes)

    @test 20 == species.vibmodes[1].Z
    @test 20 == species.vibmodes[2].Z
    @test 20 == species.vibmodes[3].Z

    @test 1918.6 == species.vibmodes[1].theta
    @test 3382.0 == species.vibmodes[2].theta
    @test 959.0 == species.vibmodes[3].theta
    
    @test 1 == species.vibmodes[1].degen
    @test 1 == species.vibmodes[2].degen
    @test 2 == species.vibmodes[3].degen
end