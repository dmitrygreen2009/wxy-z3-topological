# Exact bond-two fixtures, no DMRG/VUMPS or production checkpoint remeasurement.
using Test
include("../src/model.jl")
function validate48_tiled()
r=JSON3.read(read("results/physics48_tiled_trial.json",String),Dict{String,Any});lat=cylinder(r["family"],r["length"],r["width"];ordering="star")
s=siteinds("S=1/2",lat.n;conserve_qns=true);ip=invperm(lat.order);os=OpSum()
for (v,i,e) in r["retained_star_leg_terms"],a=1:3
 m=ip[3(v-1)+a];g=ip[3lat.nv+e];os += -W[a,i],"S-",m,"S+",g;os += -conj(W[a,i]),"S+",m,"S-",g
end
Hret=MPO(os,s);Hfull=MPO(lat.os,s);records=[]
@testset "Tiled 2D trial: independently known retained/cut Bell pair" begin
 for (name,term,expected) in [("retained",r["retained_star_leg_terms"][1],-1/sqrt(3)),("cut",r["removed_star_leg_terms"][1],0.0)]
  v,i,e=term;m=ip[3(v-1)+1];g=ip[3lat.nv+e];pattern=fill("Dn",lat.n);pattern[m]="Up";psi=MPS(s,pattern)
  rotation=exp((pi/4)*(op("S-",s[m])*op("S+",s[g])-op("S+",s[m])*op("S-",s[g])))
  psi=apply(rotation,psi;cutoff=1e-14);@test abs(norm(psi)-1)<1e-12
  full=real(inner(psi',Hfull,psi));ret=real(inner(psi',Hret,psi));@test full≈-1/sqrt(3) atol=1e-12;@test ret≈expected atol=1e-12
  push!(records,Dict("fixture"=>name,"star_leg"=>term,"full_H_energy"=>full,"retained_H_energy"=>ret,"expected_retained_energy"=>expected,"bond_dimension"=>maxlinkdim(psi)))
 end
end
atomic_json("results/physics48_tiled_trial_validation.json",Dict("records"=>records,"optimization_performed"=>false,"fixture_scope"=>"Exact one-excitation Bell states on the original shared-edge graph, testing physical-site mapping and zero expectation across a cut"))

end
validate48_tiled()
