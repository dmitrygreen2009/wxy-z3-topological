# ARCHIVAL ONLY: no distinguished microscopic diagnostic role established.
# Excluded from core phase interpretation and convergence judgments.
using Test
include("../scripts/observables.jl")
@testset "Matter triplet correlator" begin
    sites=siteinds("S=1/2",6;conserve_qns=true)
    p=MPS(sites,["Up","Up","Up","Dn","Dn","Dn"])
    q=MPS(sites,["Dn","Dn","Dn","Up","Up","Up"])
    psi=p+q;normalize!(psi)
    @test triplet_correlation(psi,[1,2,3],[4,5,6])≈0.5 atol=1e-13
    @test triplet_correlation(p,[1,2,3],[4,5,6])≈0.0 atol=1e-13
end
