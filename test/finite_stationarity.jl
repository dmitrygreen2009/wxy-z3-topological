using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/finite_stationarity.jl")
@testset "Projected stationarity of exactly checkable two-spin states" begin
 s=siteinds("S=1/2",2);os=OpSum();os+=-1.,"S+",1,"S-",2;os+=-1.,"S-",1,"S+",2;H=MPO(os,s)
 p=MPS(s,["Up","Dn"]);r=finite_bond_stationarity(H,p,1)
 @test r["projected_rayleigh_energy_real"]≈0 atol=1e-12
 @test r["projected_residual"]≈1 atol=1e-12
 v=ITensor(ComplexF64,s...);v[s[1]=>1,s[2]=>2]=1/sqrt(2);v[s[1]=>2,s[2]=>1]=1/sqrt(2)
 bell=MPS(v,s);r=finite_bond_stationarity(H,bell,1)
 @test r["projected_rayleigh_energy_real"]≈-1 atol=1e-12
 @test r["projected_residual"]<1e-12
 @test norm(bell)≈1 atol=1e-12
end
