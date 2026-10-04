isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/infinite_model.jl")
using ITensorInfiniteMPS,Test
@testset "Independent coherent-spin exchange expectation across cell periods" begin
 for family in ["zigzag","armchair"],slices in [2,3,6]
  w=1;n=9w*slices;s=infsiteinds("S=1/2",n;conserve_qns=false,initstate=j->"Dn")
  psi=InfMPS(s,j->"Dn")
  for j=1:n
   u=ITensor(ComplexF64[1 1;-1 1]/sqrt(2),prime(s[j]),dag(s[j]))
   psi.AL[j]=noprime(u*psi.AL[j]);psi.AR[j]=noprime(u*psi.AR[j])
  end
  H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering="star",cell_slices=slices),siteinds(only,psi))
  @test real(sum(expect(psi,H)))≈-sqrt(3)*w*slices atol=1e-10
  @test length(psi.AL)==n
 end
end
