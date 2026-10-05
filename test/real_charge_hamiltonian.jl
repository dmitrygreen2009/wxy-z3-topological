using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:matter_charge_basis) || include("../src/matter_charge_basis.jl")
include("../src/real_charge_hamiltonian.jl")
@testset "Exact real winding-basis microscopic MPO" begin
    Random.seed!(7257)
    for family in ["zigzag","armchair"]
        lat=cylinder(family,2,1;ordering="star")
        defs=JSON3.read(read("geometry/cgs_cycles/$(family)_L2_w1_star.json",String),Dict{String,Any})
        weights=winding_site_weights(lat,Int.(defs["cycles"][1]["edge_exponents"]))
        sites=winding_qn_sites(weights)
        original=MPO(charge_basis_opsum(lat),sites)
        realH=MPO(real_charge_basis_opsum(lat),sites)
        for q=0:2
            psi=random_mps(ComplexF64,sites,winding_initial_state(weights,2,q);linkdims=3)
            a=apply(original,psi;cutoff=1e-14,maxdim=256)
            b=apply(realH,psi;cutoff=1e-14,maxdim=256)
            @test norm(a-b)<1e-10
            @test inner(psi',realH,psi)≈inner(psi',original,psi) atol=1e-11
            @test flux(a)==flux(b)
        end
    end
end
