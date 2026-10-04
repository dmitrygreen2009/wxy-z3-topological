using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
@testset "Independent exact Schmidt spectra" begin
    s=siteinds("S=1/2",4)
    p=MPS(s,["Up","Dn","Up","Dn"])
    @test entropy_at(p,2)[1]≈0 atol=1e-14
    # Unequal Schmidt weights span full physical sites, independently of WXY.
    tensor=ITensor(s...)
    tensor[s[1]=>1,s[2]=>1,s[3]=>1,s[4]=>1]=sqrt(0.3)
    tensor[s[1]=>2,s[2]=>2,s[3]=>2,s[4]=>2]=sqrt(0.7)
    p=MPS(tensor,s;cutoff=0)
    S,probabilities=entropy_at(p,2)
    @test probabilities≈[0.7,0.3] atol=1e-14
    @test S≈-0.3log(0.3)-0.7log(0.7) atol=1e-14
    tensor=ITensor(s...)
    tensor[s[1]=>1,s[2]=>1,s[3]=>1,s[4]=>1]=1/sqrt(2)
    tensor[s[1]=>2,s[2]=>2,s[3]=>2,s[4]=>2]=1/sqrt(2)
    p=MPS(tensor,s;cutoff=0)
    @test entropy_at(p,2)[1]≈log(2) atol=1e-14
end
