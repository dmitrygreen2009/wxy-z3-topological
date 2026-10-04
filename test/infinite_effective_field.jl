using Test,ITensorInfiniteMPS
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/infinite_model.jl")
@testset "Official nonlocal VUMPS center operator matches microscopic local fields" begin
 for family in ["zigzag","armchair"],ordering in ["star","matter_first"]
  w=2;slices=2;n=9w*slices
  geometry=JSON3.read(read("geometry/infinite_$(family)_w$(w)_$(ordering).json",String),Dict{String,Any})
  alpha=[complex(cos((.7+.2sin(j))/2)) for j=1:n]
  beta=[cis(.37j+.013j^2)*sin((.7+.2sin(j))/2) for j=1:n]
  sites=infsiteinds("S=1/2",n;conserve_qns=false,initstate=j->"Up")
  psi=InfMPS(sites,j->"Up")
  for j=1:n
   U=ITensor(ComplexF64[alpha[j] -conj(beta[j]);beta[j] conj(alpha[j])],prime(sites[j]),dag(sites[j]))
   psi.AL[j]=noprime(U*psi.AL[j]);psi.AR[j]=noprime(U*psi.AR[j])
  end
  H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering=ordering),siteinds(only,psi))
  L,_=ITensorInfiniteMPS.left_environment(H,psi;tol=1e-12)
  R,_=ITensorInfiniteMPS.right_environment(H,psi;tol=1e-12)
  for j in [1,6,9,18,35,36]
   expected=zeros(ComplexF64,2,2)
   for vertex in geometry["vertices"],(a,m) in enumerate(vertex["matter_sites"]),edge in vertex["incident_gauge_sites"]
    g=edge["site_in_cell"];i=edge["leg"]
    if j==m
     c=-W[a,i]*conj(alpha[g])*beta[g];expected[2,1]+=c;expected[1,2]+=conj(c)
    elseif j==g
     c=-W[a,i]*conj(beta[m])*alpha[m];expected[1,2]+=c;expected[2,1]+=conj(c)
    end
   end
   M=ITensorInfiniteMPS.Hᴬᶜ(H,L,R,psi,j);v=psi.AL[j]*psi.C[j];ii=inds(v);dims=Tuple(dim.(ii));matrix=zeros(ComplexF64,2,2)
   for k=1:2
    data=zeros(ComplexF64,2);data[k]=1;basis=ITensor(reshape(data,dims),ii...)
    matrix[:,k]=vec(Array(M(basis),ii...))
   end
   matrix-=tr(matrix)/2*I
   @test matrix≈expected atol=1e-10
   @test norm(matrix-matrix')<1e-10
  end
 end
end
