@testset "Monte Carlo Tests" begin
    @test size(mc_supersonic(10)) == (10, 7)
end;