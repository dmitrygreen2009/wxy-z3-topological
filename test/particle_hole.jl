using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:particle_hole_mps) || include("../src/particle_hole.jl")
@testset "Exact particle-hole state mapping retains QNs and entropy" begin
    Random.seed!(7119)
    s=siteinds("S=1/2",6;conserve_qns=true)
    psi=random_mps(ComplexF64,s,["Up","Up","Dn","Dn","Dn","Dn"];linkdims=8)
    H=microscopic_mpo([[1,2,3]],s);mirror=particle_hole_mps(psi)
    @test hasqns(mirror) && val(flux(mirror),"Sz")==-val(flux(psi),"Sz")
    @test norm(mirror)≈norm(psi) atol=1e-12
    @test inner(mirror',H,mirror)≈inner(psi',H,psi) atol=1e-12
    @test entropy_at(deepcopy(mirror),3)[1]≈entropy_at(deepcopy(psi),3)[1] atol=1e-12
    @test abs(inner(psi,particle_hole_mps(mirror))-1)<1e-12
end
