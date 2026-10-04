using Test,ITensorInfiniteMPS,KrylovKit
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/transfer_analysis.jl")
@testset "Transfer length conversion is invariant under regrouping periods" begin
 p=.8
 for slices in [2,3,6]
  s=infsiteinds("S=1/2",slices;conserve_qns=false,initstate=j->"Up")
  bonds=[Index(2,"Link,l=$(j),c=1") for j=1:slices];tensors=ITensor[]
  for j=1:slices
   previous=j==1 ? ITensorInfiniteMPS.translatecelltags(bonds[end],-1) : bonds[j-1]
   bond=bonds[j];a=ITensor(ComplexF64,previous,s[j],bond)
   a[previous=>1,s[j]=>1,bond=>1]=sqrt(p)
   a[previous=>2,s[j]=>1,bond=>2]=sqrt(1-p)
   a[previous=>1,s[j]=>2,bond=>2]=sqrt(p)
   a[previous=>2,s[j]=>2,bond=>1]=sqrt(1-p)
   push!(tensors,a)
  end
  spectrum=transfer_spectrum(InfiniteMPS(tensors);tol=1e-12)
  magnitudes=sort(abs.(spectrum["values"]);rev=true)
  @test magnitudes[1:2]≈[1.,.8^slices] atol=1e-10
  @test slices*spectrum["xi_cells"]≈-1/log(.8) atol=1e-9
  @test maximum(spectrum["residuals"])<1e-10
 end
end
