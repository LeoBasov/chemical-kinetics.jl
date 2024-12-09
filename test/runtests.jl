using ChemicalKinetics
using Test

@testset "Solver.jl" begin
    clear!()

    add_species!("../data/CH4.json")
    set_nrho!(1e22)
    set_T!(300)
    set_frac!("CH4", 1.0)

    state = ChemicalKinetics.state

    nu = ChemicalKinetics.calc_coll_freq(state.species["CH4"], state.nrho, state.T)
    ekin_rot = ChemicalKinetics.calc_ekin_rot(state.T, state.mole_fractions, state.species)
    T = ChemicalKinetics.calc_Tkin_rtot(ekin_rot, state.mole_fractions, state.species)

    @test 3.6990640842564846e6 == nu
    @test 1.2425841000000001e-20 == ekin_rot
    @test 300 == T
end

@testset "ChemicalKinetics.jl" begin
    clear!()

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