using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:cgs_cycle_gates) || include("../src/cgs.jl")
isdefined(Main,:matter_charge_basis) || include("../src/matter_charge_basis.jl")
@testset "Direct winding QN is an exact microscopic basis change" begin
    Random.seed!(7253)
    for family in ["zigzag","armchair"]
        lat=cylinder(family,2,1;ordering="star")
        table=JSON3.read(read("geometry/cgs_cycles/$(family)_L2_w1_star.json",String),Dict{String,Any})
        p=Int.(table["cycles"][1]["edge_exponents"])
        weights=winding_site_weights(lat,p);sites=winding_qn_sites(weights)
        H=MPO(charge_basis_opsum(lat),sites)
        for charge in 0:2
            state=winding_initial_state(weights,2,charge)
            phi=random_mps(ComplexF64,sites,state;linkdims=3)
            @test val(flux(phi),"Sz")==4-lat.n
            @test mod(val(flux(phi),"Winding"),3)==charge
            physical=physical_from_charge_basis(phi,lat;maxdim=256)
            ps=siteinds(physical);Hphysical=MPO(lat.os,ps)
            @test norm(physical)≈1 atol=1e-11
            @test inner(phi',H,phi)≈inner(physical',Hphysical,physical) atol=1e-10
            rotated=apply(cgs_cycle_gates(lat,ps,p),physical;cutoff=1e-14,maxdim=256)
            @test inner(physical,rotated)≈cis(2pi*charge/3) atol=1e-11
            @test sum(expect(physical,"Sz"))≈2-lat.n/2 atol=1e-11
            Sbasis,_=entropy_at(phi,lat.cut);Sphysical,_=entropy_at(physical,lat.cut)
            @test Sbasis≈Sphysical atol=1e-10
        end
    end
end
