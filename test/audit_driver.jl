using Test
include("../scripts/infinite.jl")
@testset "Audited driver uses official VUMPS iteration" begin
    Random.seed!(8130)
    initstate(j)=isodd(j) ? "Up" : "Dn"
    sites=infsiteinds("S=1/2",18;conserve_qns=false,initstate)
    p=InfMPS(sites,initstate)
    for j=1:18
        theta=0.2+rand()*1.1;phase=cis(2pi*rand())
        rotation=ComplexF64[cos(theta) -conj(phase)*sin(theta);phase*sin(theta) cos(theta)]
        u=ITensor(rotation,sites[j]',dag(sites[j]))
        p.AL[j]=noprime(u*p.AL[j]);p.AR[j]=noprime(u*p.AR[j])
    end
    initial_h=InfiniteSum{MPO}(infinite_opsum("zigzag",1),sites)
    p=subspace_expansion(p,initial_h;cutoff=1e-10,maxdim=2)
    # Serialization copies all library state without its broken deepcopy method.
    buffer=IOBuffer();serialize(buffer,p);payload=take!(buffer)
    fresh()=deserialize(IOBuffer(payload))
    s=siteinds(only,p);H=InfiniteSum{MPO}(infinite_opsum("zigzag",1),s)
    Random.seed!(8131)
    reference=tdvp(H,fresh();time_step=-Inf,maxiter=1,tol=1e-7,solver_tol=x->1e-10,outputlevel=0,multisite_update_alg="parallel")
    Random.seed!(8131)
    candidate=audited_vumps(H,fresh();family="zigzag",w=1,cap=2,maxiter=1,
        tol=1e-7,solver_tol=x->1e-10,tag="_driver_parallel_local_audit",seed=8131,update_algorithm="parallel")
    @test real(sum(expect(candidate,H)))≈real(sum(expect(reference,H))) atol=1e-11
    @test maximum(abs(norm(candidate.C[j])-norm(reference.C[j])) for j=1:18)<1e-12
    raw=JSON3.read(last(readlines("results/infinite_zigzag_w1_chi2_driver_parallel_local_audit_local_eigensolves.jsonl")),Dict{String,Any})["solves"]
    @test length(raw)==36
    @test all(r["requested_tolerance"]==1e-10 for r in raw)
    @test all(haskey(r,"true_eigenpair_residual") && haskey(r,"converged_eigenpairs") for r in raw)
end
