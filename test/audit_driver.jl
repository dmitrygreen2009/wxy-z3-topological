using Test
include("../scripts/infinite.jl")
@testset "Audited driver uses official VUMPS iteration" begin
    p=deserialize("results/infinite_zigzag_w1_chi2.jls")
    s=siteinds(only,p);H=InfiniteSum{MPO}(infinite_opsum("zigzag",1),s)
    Random.seed!(8131)
    reference=tdvp(H,deserialize("results/infinite_zigzag_w1_chi2.jls");time_step=-Inf,maxiter=1,tol=1e-7,solver_tol=x->1e-10,outputlevel=0)
    Random.seed!(8131)
    candidate=audited_vumps(H,deserialize("results/infinite_zigzag_w1_chi2.jls");family="zigzag",w=1,cap=2,maxiter=1,
        tol=1e-7,solver_tol=x->1e-10,tag="_driver_regression",seed=8131)
    @test real(sum(expect(candidate,H)))≈real(sum(expect(reference,H))) atol=1e-11
    @test maximum(abs(norm(candidate.C[j])-norm(reference.C[j])) for j=1:18)<1e-12
end
