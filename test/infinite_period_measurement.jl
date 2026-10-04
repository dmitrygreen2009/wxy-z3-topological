include("../scripts/infinite.jl")
using Test
@testset "Six-slice full-state measurement units and cut accounting" begin
 slices=6;w=1;n=9w*slices
 s=infsiteinds("S=1/2",n;conserve_qns=false,initstate=j->"Dn")
 psi=InfMPS(s,j->"Dn")
 for j=1:n
  u=ITensor(ComplexF64[1 1;-1 1]/sqrt(2),prime(s[j]),dag(s[j]))
  psi.AL[j]=noprime(u*psi.AL[j]);psi.AR[j]=noprime(u*psi.AR[j])
 end
 H=InfiniteSum{MPO}(infinite_opsum("zigzag",w;ordering="star",cell_slices=slices),siteinds(only,psi))
 measure_infinite(psi,H,"zigzag",w,1,0;ordering="star",tag="_cell6_product_measurement_validation")
 r=JSON3.read(read("results/infinite_zigzag_w1_chi1_cell6_product_measurement_validation.json",String),Dict{String,Any})
 @test r["validation_fixture"]==true
 @test r["cell_slices"]==6
 @test length(r["spatial_entropies_at_slice_cuts"])==6
 @test maximum(abs.(r["spatial_entropies_at_slice_cuts"]))<1e-12
 @test r["energy_per_vertex"]≈-sqrt(3)/2 atol=1e-10
 @test r["isolated_A_star_variational_upper_bound_cell"]≈-2.4030921*6 atol=1e-12
 @test r["xi_slices"]==0
 @test occursin("cell6",r["checkpoint_file"])
end
