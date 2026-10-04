using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
@testset "Finite checkpoint observer preserves DMRG and resumes" begin
    Random.seed!(8131)
    sites=siteinds("S=1/2",6;conserve_qns=true)
    H=microscopic_mpo([[1,2,3]],sites)
    initial=random_mps(sites,["Up","Up","Up","Dn","Dn","Dn"];linkdims=8)
    reference,_=dmrg(H,deepcopy(initial);nsweeps=4,maxdim=8,cutoff=1e-13,noise=0,outputlevel=0)
    lat=(n=6,cut=3,circumference=0.0)
    audit=run_provenance(;seed=8131,solver="checkpoint regression DMRG",settings=Dict("sweeps"=>4,"cutoff"=>1e-13),
        initialization="random fixed N_up=3",conserved_quantum_numbers=["N_up=3"])
    observer=finite_checkpoint_observer(lat;family="star_benchmark",L=1,w=0,cap=8,seed=8131,
        stage="regression",ordering="matter_first",audit,cutoff=1e-13,noise=0)
    candidate,psi=dmrg(H,deepcopy(initial);nsweeps=4,maxdim=8,cutoff=1e-13,noise=0,outputlevel=0,observer)
    @test candidate≈reference atol=1e-12
    @test candidate≈-2.403092127540 atol=1e-10
    path="results/checkpoints/star_benchmark_w0_L1_N6_chi8_Nup3_seed8131_v3_matter_first_regression_latest.jls"
    loaded,metadata=load_valid_checkpoint(path)
    @test metadata["sweep"]==4
    @test entropy_at(loaded,3)[1]≈entropy_at(psi,3)[1] atol=1e-12
    @test_throws AssertionError save_checkpoint(path*".invalid",psi,merge(copy(metadata),Dict("nup"=>2)))
    @test !isfile(path*".invalid")
    @test resolve_checkpoint(path*".previous")==path*".previous"
    @test abs(inner(load_state(path),psi)-1)<1e-12
    mktempdir() do directory
        corrupt=joinpath(directory,"corrupt.jls")
        bytes=read(path);bytes[end]=xor(bytes[end],0x01);write(corrupt,bytes)
        cp(replace(path,".jls"=>".json"),replace(corrupt,".jls"=>".json"))
        @test_throws AssertionError load_valid_checkpoint(corrupt)
    end
    result=Dict("family"=>"star_benchmark","length"=>1,"width"=>0,"spins"=>6,
        "physical_circumference"=>0.0,"seed"=>8131,"nup"=>3,"ordering"=>"matter_first")
    first_saved=completed_finite_checkpoint(psi,result,8,"immutable_regression")
    second_saved=completed_finite_checkpoint(psi,result,8,"immutable_regression")
    @test first_saved["checkpoint_file"]!=second_saved["checkpoint_file"]
    @test isfile(first_saved["checkpoint_file"]) && isfile(second_saved["checkpoint_file"])
    @test abs(inner(load_state(first_saved["checkpoint_file"]),psi)-1)<1e-12
    resumed,_=dmrg(H,loaded;nsweeps=2,maxdim=8,cutoff=1e-13,noise=0,outputlevel=0)
    @test resumed≈candidate atol=1e-12
end
