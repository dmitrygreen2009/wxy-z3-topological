using Test,LinearAlgebra,Random,KrylovKit
include("../src/local_eigenpair_audit.jl")
@testset "Local eigenpair evidence detects achieved and missed tolerances" begin
    matrix=Diagonal([-2.,1.]);M(v)=matrix*v
    vals,vecs,info=KrylovKit.eigsolve(M,[1.,1.],1,:SR;ishermitian=true,tol=1e-12,krylovdim=2)
    r=local_eigenpair_record(M,vals[1],vecs[1],info,1e-12)
    @test r["local_criterion_passed"]
    @test r["true_eigenpair_residual"]<1e-12
    @test r["converged_eigenpairs"]>=1
    matrix2=Diagonal(collect(1.:30.));F(v)=matrix2*v
    vals,vecs,info=KrylovKit.eigsolve(F,ones(30),1,:SR;ishermitian=true,tol=1e-14,krylovdim=2,maxiter=1)
    r=local_eigenpair_record(F,vals[1],vecs[1],info,1e-14)
    @test !r["local_criterion_passed"]
    @test r["converged_eigenpairs"]==0
    @test r["true_eigenpair_residual"]>1e-14
end
