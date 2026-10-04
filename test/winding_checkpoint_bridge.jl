using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:cgs_cycle_gates) || include("../src/cgs.jl")
isdefined(Main,:bridge_physical_checkpoint) || include("../src/winding_checkpoint_bridge.jl")
@testset "Existing physical MPS to exact winding-QN checkpoint bridge" begin
    Random.seed!(7255)
    for family in ["zigzag","armchair"]
        lat=cylinder(family,2,1;ordering="star")
        definitions=JSON3.read(read("geometry/cgs_cycles/$(family)_L2_w1_star.json",String),Dict{String,Any})
        p=Int.(definitions["cycles"][1]["edge_exponents"])
        s=siteinds("S=1/2",lat.n;conserve_qns=true)
        psi=random_mps(ComplexF64,s,[j<=2 ? "Up" : "Dn" for j=1:lat.n];linkdims=3)
        for charge in 0:2
            reference,record=cgs_project(psi,lat,p;charge,maxdim=256,cutoff=1e-14)
            converted,bridge=bridge_physical_checkpoint(psi,lat,p;charge,maxdim=256,cutoff=1e-14)
            @test val(flux(converted),"Sz")==4-lat.n
            @test mod(val(flux(converted),"Winding"),3)==charge
            @test bridge["operator_maxlinkdim"]<=12
            @test bridge["projection_weight"]≈record["projection_weight"] atol=1e-11
            physical=physical_from_charge_basis(converted,lat;maxdim=256)
            dense_reference=dense(reference)
            dense_reference=replace_siteinds(dense_reference,siteinds(physical))
            @test abs(inner(dense_reference,physical)-1)<1e-10
            Hrot=MPO(charge_basis_opsum(lat),siteinds(converted))
            Hold=MPO(lat.os,siteinds(physical))
            @test inner(converted',Hrot,converted)≈inner(physical',Hold,physical) atol=1e-10
        end
    end
end
