using Test,ITensors,ITensorMPS,Random
include("../src/infinite_observables.jl")
@testset "Independent local-loop operator checks" begin
    s=siteinds("S=1/2",9)
    psi=MPS(s,["Up","Up","Up","Dn","Up","Dn","Dn","Dn","Dn"])
    U=infinite_twocycle_mpo(s)
    localpsi=MPS(s[1:8],["Up","Up","Up","Dn","Up","Dn","Dn","Dn"])
    @test inner(localpsi',U,localpsi)≈cis(4pi/3) atol=1e-12
    transformed=apply(U,apply(U,apply(U,localpsi;cutoff=1e-14);cutoff=1e-14);cutoff=1e-14)
    @test abs(inner(localpsi,transformed)-1)<1e-12
    Random.seed!(7117)
    random_state=random_mps(ComplexF64,s[1:8];linkdims=8)
    rotated=apply(U,apply(U,apply(U,random_state;cutoff=1e-14);cutoff=1e-14);cutoff=1e-14)
    @test abs(inner(random_state,rotated)-1)<1e-11
end
using ITensorInfiniteMPS
@testset "Legacy local-loop entry point translates scalar indices" begin
    initstate(j)=mod1(j,18)==4 ? "Up" : "Dn"
    sites=infsiteinds("S=1/2",18;initstate)
    psi=InfMPS(sites,initstate)
    for offset in [0,9,18]
        expected=offset==9 ? 1.0 : cis(2pi/3)
        @test expect(psi,infinite_twocycle_mpo(sites,offset))≈expected atol=1e-11
    end
end
