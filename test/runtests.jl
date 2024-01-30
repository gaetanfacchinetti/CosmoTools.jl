using CosmoTools
using Test

@testset "CosmoTools.jl" begin
    @test all(CosmoTools.Ω.(CosmoTools.Matter, 10. .^range(-2, 6, 100)) .<= 1.0)
    @test all(CosmoTools.Ω.(CosmoTools.Radiation, 10. .^range(-2, 6, 100)) .<= 1.0)
    @test all(CosmoTools.Ω.(CosmoTools.DarkEnergy, 10. .^range(-2, 6, 100)) .<= 1.0)
    @test all(CosmoTools.transfer_function.(10 .^range(-4, 6, 100)) .>= 0)
    @test all(CosmoTools.matter_power_spectrum.(10 .^range(-4, 6, 100)) .>= 0)
    @test all(CosmoTools.matter_power_spectrum.(10 .^range(-4, 6, 100)) .>= 0)
end
