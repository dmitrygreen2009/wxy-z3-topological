using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
@testset "Individual-site star ordering realizes the declared spatial partition" begin
    for family in ["zigzag","armchair"], L in [2,3,4,5], w=1:3
        lat=cylinder(family,L,w;ordering="star")
        gaugekeys=fill(NaN,lat.n-3lat.nv)
        for (v,(t,j,s)) in enumerate(lat.verts),i=1:3
            at=s==0 ? t : i==1 ? t : i==2 ? t+1 : family=="armchair" ? t-1 : t
            bt=i==1 ? at : i==2 ? at-1 : family=="armchair" ? at+1 : at
            gaugekeys[lat.legs[v][i]]=(at+bt)/2+0.1
        end
        keys=vcat([t+0.01a for (t,j,s) in lat.verts for a=1:3],gaugekeys)
        expected=findall(k->k<L/2,keys)
        @test sort(lat.order[1:lat.cut])==expected
        if iseven(L)
            # Existing even-length production ordering remains exactly the same.
            axial=cylinder(family,L,w;ordering="axial")
            @test sort(axial.order[1:axial.cut])==expected
        end
    end
end
