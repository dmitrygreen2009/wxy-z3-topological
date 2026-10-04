using Test,ITensors,ITensorMPS,ITensorInfiniteMPS
isdefined(Main,:infinite_triplet_pair_mpo) || include("../src/infinite_observables.jl")
@testset "Charge-three correlations against exact product expectations" begin
    s=siteinds("S=1/2",12)
    plus=MPS([ITensor([1.,1.]/sqrt(2),j) for j in s])
    @test inner(plus',infinite_triplet_pair_mpo(s,1,10),plus)≈1/64 atol=1e-12
    down=MPS(s,fill("Dn",12))
    @test inner(down',infinite_triplet_pair_mpo(s,1,10),down)≈0 atol=1e-12
    initstate(j)="Dn"
    si=infsiteinds("S=1/2",9;initstate)
    product=InfMPS(si,initstate)
    @test expect(product,infinite_triplet_pair_mpo(si,1,19))≈0 atol=1e-12
    for j=1:9
        h=op("H",si[j])
        product.AL[j]=apply(h,product.AL[j])
        product.AR[j]=apply(h,product.AR[j])
    end
    @test expect(product,infinite_triplet_pair_mpo(si,1,19))≈1/64 atol=1e-12
    # H|Dn>=(|Up>-|Dn>)/sqrt(2), so each S+ averages to -1/2.
    @test expect(product,infinite_triplet_mpo(si,1))≈-1/8 atol=1e-12
end
