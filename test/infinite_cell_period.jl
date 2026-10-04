isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/infinite_model.jl")
using Test
@testset "Individual-site microscopic infinite-cell periods" begin
 for family in ["zigzag","armchair"], w in [1,2,3], order in ["star","matter_first"], slices in [2,3,6]
  n=9w*slices
  pairs=infinite_pairs(family,w;ordering=order,cell_slices=slices)
  @test length(pairs)==18w*slices
  @test all(p->1<=min(p.m,p.g)<=n && p.m!=p.g,pairs)
  counts=zeros(Int,n)
  for p in pairs
   counts[mod1(p.m,n)]+=1;counts[mod1(p.g,n)]+=1
  end
  # Each matter spin couples to three physical gauges; each gauge to six matter spins.
  @test count(==(3),counts)==6w*slices
  @test count(==(6),counts)==3w*slices
  @test length(unique((p.m,p.g,p.a,p.i) for p in pairs))==length(pairs)
  @test length(infinite_opsum(family,w;ordering=order,cell_slices=slices))==2length(pairs)+n
 end
end
