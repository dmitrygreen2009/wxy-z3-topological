using Test
include("../src/resume_controls.jl")
@testset "Resume retains recorded local solver accuracy" begin
    @test Meta.parseall(read(joinpath(@__DIR__,"../scripts/resume_checkpoint.jl"),String)) isa Expr
    strict=Dict("krylovdim"=>32,"eigsolve_maxiter"=>10,"eigsolve_tol"=>1e-14,"eigsolve_verbosity"=>1)
    c=finite_resume_controls(strict,"continuation_pass1")
    @test c==(eigsolve_krylovdim=32,eigsolve_maxiter=10,eigsolve_tol=1e-14,eigsolve_verbosity=1)
    staged=Dict("krylovdim"=>8,"refinement_krylovdim"=>12)
    @test finite_resume_controls(staged,"cap256").eigsolve_krylovdim==8
    @test finite_resume_controls(staged,"refinement").eigsolve_krylovdim==12
    @test finite_resume_controls(Dict(),"legacy").eigsolve_tol==1e-14
end
