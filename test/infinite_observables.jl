using Test,ITensors,ITensorMPS,Random
include("../src/infinite_observables.jl")
@testset "Independent local-loop and matter-triplet operator checks" begin
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
    # Independent tensor-product expectation for a three-spin coherent state.
    triplet_sites=siteinds("S=1/2",3)
    product=MPS([ITensor(ComplexF64[1,1]/sqrt(2),j) for j in triplet_sites])
    @test inner(product',infinite_triplet_mpo(triplet_sites,1),product)≈1/8 atol=1e-12
end
