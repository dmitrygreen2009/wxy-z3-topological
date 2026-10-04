using Test,ITensorInfiniteMPS
isdefined(Main,:cylinder) || include("../src/model.jl")
@testset "Infinite bond dimensions include the larger wrap bond" begin
 s=infsiteinds("S=1/2",3;conserve_qns=false,initstate=j->"Up")
 links=[Index(d,"Link,l=$(j),c=1") for (j,d) in enumerate([2,2,5])]
 tensors=ITensor[]
 for j=1:3
  previous=j==1 ? ITensorInfiniteMPS.translatecelltags(links[end],-1) : links[j-1]
  push!(tensors,ITensor(ComplexF64,previous,s[j],links[j]))
 end
 A=InfiniteMPS(tensors);details=mps_bond_dimension_details(A)
 @test prod(dim.(dag(input_inds(TransferMatrix(A)))))==25
 @test maximum_bond_dimension(A)==5
 @test details["bond_dimensions_by_tensor_representation"]["MPS"]==[2,2,5]
 @test details["periodic_wrap_bond_included"]
 C=InfiniteMPS([ITensor() for _=1:3]);canonical=InfiniteCanonicalMPS(A,C,A)
 @test maximum_bond_dimension(canonical)==5
 @test mps_bond_dimension_details(canonical)["bond_dimensions_by_tensor_representation"]["AR"]==[2,2,5]
 finite=MPS(siteinds("S=1/2",3),["Up","Dn","Up"])
 @test maximum_bond_dimension(finite)==maxlinkdim(finite)==1
 @test !mps_bond_dimension_details(finite)["periodic_wrap_bond_included"]
end
