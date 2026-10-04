using Test,ITensorInfiniteMPS
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/cgs.jl")
include("../src/cgs_operator_mpo.jl")
@testset "Factorized microscopic circumference MPO and translated infinite support" begin
    Random.seed!(7221)
    for family in ["zigzag","armchair"],w in [1,2]
        lat=cylinder(family,4,w;ordering="star");s=siteinds("S=1/2",lat.n)
        leading=lat.n-36w;bulk=s[leading+1:leading+18w]
        spec=cgs_circumference_spec(family,w);U=cgs_product_mpo(bulk,spec)
        psi=random_mps(ComplexF64,bulk[spec.first_site:spec.last_site];linkdims=3)
        table=JSON3.read(read("geometry/cgs_cycles/$(family)_L4_w$(w)_star.json",String),Dict{String,Any})
        p=table["cycles"][1]["edge_exponents"]
        reference=apply(cgs_cycle_gates(lat,s,p),psi;cutoff=1e-14,maxdim=128)
        measured=apply(U,psi;cutoff=1e-14,maxdim=128)
        @test norm(measured)≈1 atol=1e-11
        @test abs(inner(reference,measured)-1)<1e-11
        @test maxlinkdim(U)<=4
        n=18w;initstate(j)=mod1(j,n)==4 ? "Up" : "Dn"
        infinite_sites=infsiteinds("S=1/2",n;initstate)
        product=InfMPS(infinite_sites,initstate)
        for offset in [0,9w,18w]
            operator=cgs_product_mpo(infinite_sites,cgs_circumference_spec(family,w;offset))
            expected=family=="zigzag" && offset==9w ? 1.0 : cis(2pi/3)
            @test expect(product,operator)≈expected atol=1e-11
        end
    end
end
