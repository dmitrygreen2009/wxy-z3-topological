using Test,ITensorInfiniteMPS,KrylovKit
isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:infinite_number_background) || include("../src/infinite_number_background.jl")
@testset "Shifted zero flux distinguishes physical fillings" begin
    for (n,nup,label) in [(18,9,"Sz0"),(36,16,"rho4of9"),(36,20,"rho5of9")]
        initstate(j)=fld(j*nup,n)>fld((j-1)*nup,n) ? "Up" : "Dn"
        s=infsiteinds("S=1/2",n;conserve_qns=true,initstate)
        psi=InfMPS(s,initstate);background=infinite_number_background(psi)
        @test background["label"]==label
        @test background["target_average_nup_per_cell"]==nup
        @test background["physical_number_density"]≈nup/n atol=1e-15
        @test sum(real(expect(psi,"Sz",j))+.5 for j=1:n)≈nup atol=1e-12
        for opname in ("Sz","S+","S-")
            ordinary=siteind("S=1/2")
            @test Array(op(opname,s[1]),prime(s[1]),dag(s[1]))≈Array(op(opname,ordinary),prime(ordinary),dag(ordinary)) atol=1e-15
        end
    end
end
