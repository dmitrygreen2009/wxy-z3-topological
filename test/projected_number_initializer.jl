using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/projected_number_initializer.jl")
isdefined(Main,:cgs_product_mpo) || include("../src/cgs_operator_mpo.jl")
@testset "Unequal winding weights restrict two-site randomization" begin
    for left in 0:2,right in 0:2
        left==right && continue
        s=winding_qn_sites([left,right])
        gate=ITensorMPS.randomU(Float64,s[1],s[2])
        matrix=reshape(Array(gate,s[1]',s[2]',dag(s[1]),dag(s[2])),4,4)
        @test norm(matrix-Diagonal(diag(matrix)))<1e-14
    end
end
@testset "Physical-number initialization projects into every target winding block" begin
    for family in ["zigzag","armchair"],q in 0:2
        lat=cylinder(family,2,1;ordering="star")
        defs=JSON3.read(read("geometry/cgs_cycles/$(family)_L2_w1_star.json",String),Dict{String,Any})
        cycle=first(c for c in defs["cycles"] if c["winding_number"]==1)
        nup=family=="zigzag" ? 9 : 8
        psi,record=projected_number_initializer(lat,Int.(cycle["edge_exponents"]),nup,q;linkdims=2,maxdim=64)
        @test val(flux(psi),"Sz")==2nup-lat.n
        @test mod(val(flux(psi),"Winding"),3)==q
        @test norm(psi)≈1. atol=1e-11
        @test record["projection_weight"]>1e-12
        weights=record["onsite_winding_weights"]
        U=cgs_product_mpo(siteinds(psi),(;triplets=Tuple{Int,Int,Int}[],gauges=collect(enumerate(weights)),first_site=1,last_site=lat.n))
        @test inner(psi',U,psi)≈cis(2pi*q/3) atol=1e-10
    end
end
