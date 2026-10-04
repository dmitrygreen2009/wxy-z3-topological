using Test,ITensorInfiniteMPS
isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:cgs_cycle_gates) || include("../src/cgs.jl")
isdefined(Main,:cgs_product_mpo) || include("../src/cgs_operator_mpo.jl")
@testset "Contractible microscopic operators from explicit lifted geometry" begin
    Random.seed!(7223)
    for family in ["zigzag","armchair"],w in [1,2]
        lat=cylinder(family,4,w;ordering="star");s=siteinds("S=1/2",lat.n)
        table=JSON3.read(read("geometry/cgs_plaquettes/$(family)_L4_w$(w)_star.json",String),Dict{String,Any})
        cycle=first(r for r in table["plaquettes"] if haskey(r,"infinite_star_operator"))
        r=cycle["infinite_star_operator"]
        spec=(triplets=[Tuple(Int.(x)) for x in r["triplets"]],gauges=[Tuple(Int.(x)) for x in r["gauges"]],first_site=r["first_site"],last_site=r["last_site"])
        leading=lat.n-36w;bulk=s[leading+1:end]
        U=cgs_product_mpo(bulk,spec)
        psi=random_mps(ComplexF64,bulk[spec.first_site:spec.last_site];linkdims=2)
        reference=apply(cgs_cycle_gates(lat,s,cycle["edge_exponents"]),psi;cutoff=1e-14,maxdim=128)
        rotated=apply(U,psi;cutoff=1e-14,maxdim=128)
        @test abs(inner(reference,rotated)-1)<1e-11
        @test norm(rotated)≈1 atol=1e-11
        initstate(j)="Dn";si=infsiteinds("S=1/2",18w;initstate)
        product=InfMPS(si,initstate)
        @test expect(product,cgs_product_mpo(si,spec))≈1 atol=1e-11
    end
end
