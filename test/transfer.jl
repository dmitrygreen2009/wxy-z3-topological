using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
using ITensorInfiniteMPS, KrylovKit
include("../src/transfer_analysis.jl")
include("../src/dense_transfer_matrix.jl")
include("../src/transfer_uniqueness.jl")
@testset "Analytically known infinite transfer spectrum" begin
    # A^up=diag(sqrt(p),sqrt(1-p)); A^down swaps the two bond states.
    # The channel eigenvalues are 1, 2sqrt(p(1-p)), 0, 0.
    p=0.8;s=infsiteinds("S=1/2",1;conserve_qns=false,initstate=n->"Up")
    bond=Index(2,"Link,l=1,c=1")
    previous=ITensorInfiniteMPS.translatecelltags(bond,-1)
    tensor=ITensor(ComplexF64,previous,s[1],bond)
    tensor[previous=>1,s[1]=>1,bond=>1]=sqrt(p)
    tensor[previous=>2,s[1]=>1,bond=>2]=sqrt(1-p)
    tensor[previous=>1,s[1]=>2,bond=>2]=sqrt(p)
    tensor[previous=>2,s[1]=>2,bond=>1]=sqrt(1-p)
    A=InfiniteMPS([tensor]);T=TransferMatrix(A)
    certificate=transfer_uniqueness_certificate(A)
    @test certificate["fixed_point_uniqueness_certified"]
    @test certificate["primitive_peripheral_spectrum_certified"]
    @test !transfer_uniqueness_certificate(A;maximum_dimension=1)["fixed_point_uniqueness_certified"]
    dense_values=eigvals(dense_transfer_matrix(T))
    @test sort(abs.(dense_values);rev=true)≈[1.,.8,0.,0.] atol=1e-10
    x=random_itensor(ComplexF64,dag(input_inds(T)))
    values,vectors,info=eigsolve(T,x,4,:LM;tol=1e-12,krylovdim=8)
    magnitudes=sort(abs.(values);rev=true)
    # The population transfer has eigenvalues 1 and 0; coherence 0.8 and 0.
    @test magnitudes[1:2]≈[1.0,0.8] atol=1e-10
    @test all(magnitudes[3:end].<1e-10)
    @test -1/log(magnitudes[2]/magnitudes[1])≈-1/log(0.8) atol=1e-9
    @test maximum(norm(T(vectors[k])-values[k]*vectors[k]) for k in eachindex(values))<1e-10
end

@testset "Degenerate fixed points are not hidden by scalar Krylov" begin
    s=infsiteinds("S=1/2",1;conserve_qns=false,initstate=n->"Up")
    bond=Index(2,"Link,l=1,c=1");previous=ITensorInfiniteMPS.translatecelltags(bond,-1)
    tensor=ITensor(ComplexF64,previous,s[1],bond)
    tensor[previous=>1,s[1]=>1,bond=>1]=1
    tensor[previous=>2,s[1]=>2,bond=>2]=1
    certificate=transfer_uniqueness_certificate(InfiniteMPS([tensor]))
    @test !certificate["fixed_point_uniqueness_certified"]
    @test !certificate["primitive_peripheral_spectrum_certified"]
    @test certificate["fixed_point_multiplicity_at_tolerance"]==2
    @test sort(abs.(eigvals(dense_transfer_matrix(TransferMatrix(InfiniteMPS([tensor])))));rev=true)≈[1.,1.,0.,0.] atol=1e-11
    spectrum=transfer_spectrum(InfiniteMPS([tensor]);tol=1e-12)
    @test spectrum["fixed_point_rank_detected"]==2
    @test isinf(spectrum["xi_cells"])
    @test sort(abs.(spectrum["values"]);rev=true)[1:2]≈[1,1] atol=1e-11
end
