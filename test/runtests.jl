using Test
include("../src/model.jl")
include("../src/infinite_model.jl")
@testset "Microscopic lattice regression" begin
    @test norm(W'*W-I)<1e-14
    for family in ["zigzag","armchair"], L in [2,4], w in [1,2,3]
        lat=cylinder(family,L,w)
        @test lat.nv==2L*w
        @test all(length(l)==3 for l in lat.legs)
        @test all(length(unique(l))==3 for l in lat.legs)
        multiplicities=[count(==(e),vcat(lat.legs...)) for e=1:maximum(vcat(lat.legs...))]
        @test all(m in [1,2] for m in multiplicities)
        @test sum(multiplicities)==3lat.nv
        @test sort(lat.order)==collect(1:lat.n)
        @test 0<lat.cut<lat.n
        star=cylinder(family,L,w;ordering="star")
        @test Set(star.order[1:star.cut])==Set(lat.order[1:lat.cut])
    end
end
@testset "Infinite Hamiltonian physical incidence" begin
    for family=["zigzag","armchair"],w=1:3
        pairs=infinite_pairs(family,w);n=18w
        @test length(pairs)==36w
        @test all(1<=min(p.m,p.g)<=n for p in pairs)
        matter=sort(unique(mod1(p.m,n) for p in pairs))
        gauge=sort(unique(mod1(p.g,n) for p in pairs))
        @test length(matter)==12w
        @test length(gauge)==6w
        @test isempty(intersect(matter,gauge))
        @test all(count(p->mod1(p.m,n)==m,pairs)==3 for m in matter)
        @test all(count(p->mod1(p.g,n)==g,pairs)==6 for g in gauge)
    end
    for family=["zigzag","armchair"],w=1:3
        n=18w;permutation=zeros(Int,n)
        for t=0:1,j=0:w-1,a=1:3
            permutation[9w*t+6j+a]=9w*t+9j+a
            permutation[9w*t+6j+3+a]=9w*t+9j+5+a
        end
        for t=0:1,j=0:w-1
            permutation[9w*t+6w+3j+1]=9w*t+9j+4
            permutation[9w*t+6w+3j+2]=9w*t+9j+5
            permutation[9w*t+6w+3j+3]=9w*t+9j+9
        end
        @test sort(permutation)==collect(1:n)
        translate(i)=permutation[mod1(i,n)]+fld(i-1,n)*n
        old=Set((translate(p.m),translate(p.g),p.a,p.i) for p in infinite_pairs(family,w))
        new=Set((p.m,p.g,p.a,p.i) for p in infinite_pairs(family,w;ordering="star"))
        @test old==new
    end
    # The zigzag w=1, two-slice periodic quotient is the validated 4-star torus.
    pairs=infinite_pairs("zigzag",1)
    legs=[[7,18,8],[7,9,8],[16,9,17],[16,18,17]]
    matter=[[1,2,3],[4,5,6],[10,11,12],[13,14,15]]
    @test all(any(p->mod1(p.m,18)==matter[v][a] && mod1(p.g,18)==legs[v][i] && p.a==a && p.i==i,pairs)
        for v=1:4,a=1:3,i=1:3)
end
@testset "MPO Hermiticity" begin
    Random.seed!(7133)
    for (legs,nup) in [([[1,2,3]],3),([[1,2,3],[1,4,5]],5),([[1,2,3],[1,5,3],[4,5,6],[4,2,6]],9)]
        n=3length(legs)+maximum(vcat(legs...))
        sites=siteinds("S=1/2",n;conserve_qns=true)
        H=microscopic_mpo(legs,sites)
        state=[i<=nup ? "Up" : "Dn" for i=1:n]
        p=random_mps(ComplexF64,sites,state;linkdims=8)
        q=random_mps(ComplexF64,sites,state;linkdims=8)
        @test abs(inner(p',H,q)-conj(inner(q',H,p)))<1e-12
    end
end
