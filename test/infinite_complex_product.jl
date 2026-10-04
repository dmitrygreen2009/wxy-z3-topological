using Test,ITensorInfiniteMPS
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/infinite_model.jl")
@testset "Complex local-spin phases independently check infinite microscopic MPO" begin
 for family in ["zigzag","armchair"],slices in [2,3,6]
  w=2;n=9w*slices;part=slices==2 ? "" : "_cell$(slices)"
  geometry=JSON3.read(read("geometry/infinite_$(family)_w$(w)$(part)_star.json",String),Dict{String,Any})
  alpha=[complex(cos((.7+.2sin(j))/2)) for j=1:n]
  beta=[cis(.37j+.013j^2)*sin((.7+.2sin(j))/2) for j=1:n]
  expected=0.
  for v in geometry["vertices"],(a,m) in enumerate(v["matter_sites"]),edge in v["incident_gauge_sites"]
   g=edge["site_in_cell"];i=edge["leg"]
   expected-=2real(W[a,i]*conj(beta[m])*alpha[m]*conj(alpha[g])*beta[g])
  end
  sites=infsiteinds("S=1/2",n;conserve_qns=false,initstate=j->"Up")
  psi=InfMPS(sites,j->"Up")
  for j=1:n
   U=ITensor(ComplexF64[alpha[j] -conj(beta[j]);beta[j] conj(alpha[j])],prime(sites[j]),dag(sites[j]))
   psi.AL[j]=noprime(U*psi.AL[j]);psi.AR[j]=noprime(U*psi.AR[j])
  end
  H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering="star",cell_slices=slices),siteinds(only,psi))
  actual=sum(expect(psi,H))
  @test real(actual)≈expected atol=1e-10
  @test abs(imag(actual))<1e-10
 end
end
