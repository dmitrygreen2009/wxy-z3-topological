using Test,ITensorInfiniteMPS,KrylovKit
isdefined(Main,:cylinder) || include("../src/model.jl")
@testset "Official parallel VUMPS: exact ferromagnetic Ising ground state" begin
    Random.seed!(7232)
    s=infsiteinds("S=1/2",2;conserve_qns=false,initstate=j->"Dn")
    psi=InfMPS(s,j->"Dn")
    for j=1:2
        theta=.85+.05j;u=ITensor(ComplexF64[cos(theta) -sin(theta);sin(theta) cos(theta)],prime(s[j]),dag(s[j]))
        psi.AL[j]=noprime(u*psi.AL[j]);psi.AR[j]=noprime(u*psi.AR[j])
    end
    os=OpSum();os+=-1.,"Sz",1,"Sz",2;os+=-1.,"Sz",2,"Sz",3
    H=InfiniteSum{MPO}(os,siteinds(only,psi))
    for algorithm in ["sequential","parallel"]
        epsL=fill(1e-7,2);epsR=fill(1e-7,2)
        buffer=IOBuffer();serialize(buffer,psi);candidate=deserialize(IOBuffer(take!(buffer)))
        candidate,_=ITensorInfiniteMPS.tdvp_iteration(ITensorInfiniteMPS.vumps_solver,H,candidate;
            (ϵᴸ!)=epsL,(ϵᴿ!)=epsR,multisite_update_alg=algorithm,time_step=-Inf,solver_tol=x->1e-10,eager=true)
        @test real(sum(expect(candidate,H)))≈-.5 atol=1e-11
        @test [real(expect(candidate,"Sz",j)) for j=1:2]≈[.5,.5] atol=1e-11
        @test max(maximum(epsL),maximum(epsR))<1e-11
    end
end
