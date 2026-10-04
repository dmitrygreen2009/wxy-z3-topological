using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
using ITensorInfiniteMPS,KrylovKit
include("../src/canonical_analysis.jl")
@testset "Infinite Schmidt entropy from independently known fixed point" begin
    Random.seed!(7114);p=0.8
    s=infsiteinds("S=1/2",1;conserve_qns=false,initstate=n->"Up")
    bond=Index(2,"Link,Left,l=1,c=1");previous=ITensorInfiniteMPS.translatecelltags(bond,-1)
    tensor=ITensor(ComplexF64,previous,s[1],bond)
    tensor[previous=>1,s[1]=>1,bond=>1]=sqrt(p)
    tensor[previous=>2,s[1]=>1,bond=>2]=sqrt(1-p)
    tensor[previous=>1,s[1]=>2,bond=>2]=sqrt(p)
    tensor[previous=>2,s[1]=>2,bond=>1]=sqrt(1-p)
    state,lambda=canonicalize_left(InfiniteMPS([tensor]))
    @test infinite_schmidt_entropies(state)[1]≈-p*log(p)-(1-p)*log(1-p) atol=1e-11
    @test abs(lambda-1)<1e-11
    unnamed=InfiniteMPS([prime(removetags(tensor,"Left";tags="Link");tags="Link")])
    renamed,_=canonicalize_left(unnamed)
    @test infinite_schmidt_entropies(renamed)[1]≈-p*log(p)-(1-p)*log(1-p) atol=1e-11
    @test norm(state.AL[1]*state.C[1]-state.C[0]*state.AR[1])<1e-11
end
@testset "Conserved-spin infinite entropy via measurement-only densification" begin
    Random.seed!(7114);p=0.8
    s=infsiteinds("S=1/2",2;conserve_qns=true,initstate=n->isodd(n) ? "Up" : "Dn")
    middle=Index(QN("Sz",-1)=>1,QN("Sz",1)=>1;tags="Link,Left,l=1,c=1")
    last=Index(QN("Sz",0)=>1;tags="Link,Left,l=2,c=1")
    previous=dag(ITensorInfiniteMPS.translatecelltags(last,-1))
    t1=ITensor(ComplexF64,previous,s[1],middle)
    t1[previous=>1,s[1]=>1,middle=>1]=1
    t1[previous=>1,s[1]=>2,middle=>2]=1
    l=dag(middle);t2=ITensor(ComplexF64,l,s[2],last)
    t2[l=>1,s[2]=>2,last=>1]=sqrt(p)
    t2[l=>2,s[2]=>1,last=>1]=sqrt(1-p)
    AL=InfiniteMPS([t1,t2]);state,lambda=canonicalize_left(AL)
    @test infinite_schmidt_entropies(state)≈[-p*log(p)-(1-p)*log(1-p),0] atol=1e-11
    @test hasqns(siteind(AL,1))
    @test !hasqns(siteind(state.AL,1))
end
