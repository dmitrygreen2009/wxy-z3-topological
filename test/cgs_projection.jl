using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/cgs.jl")
@testset "Exact CGS projector and Hamiltonian commutation" begin
    lat=cylinder("zigzag",2,1;ordering="star")
    sites=siteinds("S=1/2",lat.n;conserve_qns=true)
    Random.seed!(7112)
    psi=random_mps(ComplexF64,sites,[isodd(i) ? "Up" : "Dn" for i=1:lat.n];linkdims=8)
    v=findfirst(==((0,0,0)),lat.verts);p=zeros(Int,lat.n-3lat.nv)
    p[lat.legs[v][1]]=1;p[lat.legs[v][3]]=2
    gates=cgs_cycle_gates(lat,sites,p);H=MPO(lat.os,sites)
    u=apply(gates,psi;cutoff=1e-14,maxdim=256)
    u3=apply(gates,apply(gates,u;cutoff=1e-14,maxdim=256);cutoff=1e-14,maxdim=256)
    @test abs(inner(psi,u3)-1)<1e-11
    @test inner(psi',H,psi)≈inner(u',H,u) atol=1e-11
    projected,audit=cgs_project(psi,lat,p;maxdim=256,cutoff=1e-14)
    audit["edge_exponents"]=p
    @test audit["edge_exponents"]==p
    @test audit["purity_error"]<1e-11
    for charge in [1,2]
        charged,record=cgs_project(psi,lat,p;charge,maxdim=256,cutoff=1e-14)
        @test record["purity_error"]<1e-11
        @test abs(inner(projected,charged))<1e-11
    end
    second,_=cgs_project(projected,lat,p;maxdim=256,cutoff=1e-14)
    @test abs(inner(second,projected))≈1 atol=1e-11
end
