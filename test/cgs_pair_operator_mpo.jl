using Test,ITensorInfiniteMPS
isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:cgs_cycle_gates) || include("../src/cgs.jl")
isdefined(Main,:cgs_product_mpo) || include("../src/cgs_operator_mpo.jl")
@testset "Separated winding-loop adjoint products" begin
    Random.seed!(7224)
    for family in ["zigzag","armchair"],w in [1,2]
        lat=cylinder(family,4,w;ordering="star");s=siteinds("S=1/2",lat.n)
        leading=lat.n-36w;bulk=s[leading+1:end]
        spec=cgs_disjoint_pair_spec(cgs_circumference_spec(family,w),cgs_circumference_spec(family,w;offset=18w))
        U=cgs_product_mpo(bulk,spec)
        table=JSON3.read(read("geometry/cgs_cycles/$(family)_L4_w$(w)_star.json",String),Dict{String,Any})
        p=mod.(table["cycles"][1]["edge_exponents"]-table["cycles"][3]["edge_exponents"],3)
        psi=random_mps(ComplexF64,bulk[spec.first_site:spec.last_site];linkdims=2)
        reference=apply(cgs_cycle_gates(lat,s,p),psi;cutoff=1e-14,maxdim=128)
        rotated=apply(U,psi;cutoff=1e-14,maxdim=128)
        @test abs(inner(reference,rotated)-1)<1e-11
        @test norm(rotated)≈1 atol=1e-11
        initstate(j)=mod1(j,18w)==4 ? "Up" : "Dn"
        si=infsiteinds("S=1/2",18w;initstate)
        @test expect(InfMPS(si,initstate),cgs_product_mpo(si,spec))≈1 atol=1e-11
    end
end
