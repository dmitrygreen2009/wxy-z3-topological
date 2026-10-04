using Test,ITensorInfiniteMPS
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/infinite_block_initializer.jl")
include(joinpath(pkgdir(ITensorInfiniteMPS),"examples","vumps","src","vumps_subspace_expansion.jl"))
@testset "Independent-cell initializer: exact Bell entropy and physical filling" begin
    s=siteinds("S=1/2",2;conserve_qns=true)
    phi=add(MPS(s,["Up","Dn"]),MPS(s,["Dn","Up"]);cutoff=1e-14);normalize!(phi)
    psi=infinite_cell_product(phi)
    @test hasqns(siteind(psi.AL,1))
    probabilities=abs2.(svdvals(Array(dense(psi.C[1]),inds(dense(psi.C[1]))...)))
    @test -sum(p>0 ? p*log(p) : 0. for p in probabilities)≈log(2) atol=1e-12
    @test norm(psi.C[2])≈1 atol=1e-12
    @test dim(only(linkinds(psi.AL,2=>3)))==1
    @test expect(psi,"Sz",1)≈0 atol=1e-12
    si=siteinds(only,psi);os=OpSum();os+=1.,"Sz",1,"Sz",2
    @test expect(psi,MPO(os,[si[1],si[2]]))≈-1/4 atol=1e-12
    os=OpSum();os+=1.,"Sz",1,"Sz",3
    @test expect(psi,MPO(os,[si[j] for j=1:3]))≈0 atol=1e-12
    # The official shifted/scaled QN convention enforces Nup=4 per 9 spins,
    # while the physical Sz operator retains its spin-1/2 eigenvalues.
    Random.seed!(7225);initstate(j)=j<=4 ? "Up" : "Dn"
    shifted=infsiteinds("S=1/2",9;conserve_qns=true,initstate)
    random_cell=random_mps(ComplexF64,collect(shifted[Cell(1)]),[initstate(j) for j=1:9];linkdims=4)
    @test iszero(flux(random_cell))
    tiled=infinite_cell_product(random_cell)
    @test sum(real(expect(tiled,"Sz",j))+.5 for j=1:9)≈4 atol=1e-12
    @test hasqns(siteind(tiled.AL,1))
    # Exercise official expansion across the deliberate cell boundary; exact
    # expectations alone do not check the library's QN arrow convention.
    os=OpSum()
    for j=1:2
        os+=-.5,"S+",j,"S-",j+1
        os+=-.5,"S-",j,"S+",j+1
    end
    H=InfiniteSum{MPO}(os,siteinds(only,psi))
    expanded=subspace_expansion(psi,H;cutoff=1e-12,maxdim=4)
    @test hasqns(siteind(expanded.AL,1))
    @test dim(only(linkinds(expanded.AL,2=>3)))>1
    @test all(iszero(flux(expanded.AL[j])) for j=1:2)
end
