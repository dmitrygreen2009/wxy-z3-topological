using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:cgs_cycle_gates) || include("../src/cgs.jl")
isdefined(Main,:cgs_sector_penalty) || include("../src/cgs_sector_penalty.jl")
@testset "Exact winding-sector penalty on both microscopic geometries" begin
    Random.seed!(7251)
    for family in ["zigzag","armchair"]
        lat=cylinder(family,2,1;ordering="star")
        sites=siteinds("S=1/2",lat.n;conserve_qns=true)
        table=JSON3.read(read("geometry/cgs_cycles/$(family)_L2_w1_star.json",String),Dict{String,Any})
        p=table["cycles"][1]["edge_exponents"]
        psi=random_mps(ComplexF64,sites,[j<=2 ? "Up" : "Dn" for j=1:lat.n];linkdims=2)
        H=MPO(lat.os,sites)
        for charge in 0:2
            penalty=cgs_sector_penalty(lat,sites,p;charge)
            projected,record=cgs_project(psi,lat,p;charge,maxdim=256,cutoff=1e-14)
            @test record["purity_error"]<1e-11
            @test inner(projected',penalty.U,projected)≈cis(2pi*charge/3) atol=1e-11
            @test abs(sum(inner(projected',term,projected) for term in penalty.terms))<1e-10
            ref=apply(cgs_cycle_gates(lat,sites,p),psi;cutoff=1e-14,maxdim=256)
            rotated=apply(penalty.U,psi;cutoff=1e-14,maxdim=256)
            @test norm(ref-rotated)<1e-11
            @test inner(rotated',H,rotated)≈inner(psi',H,psi) atol=1e-11
            @test penalty.strength>2penalty.hamiltonian_norm_bound
        end
    end
end
