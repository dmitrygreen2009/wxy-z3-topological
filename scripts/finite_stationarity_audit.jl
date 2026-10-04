include("cylinders.jl")
include("../src/finite_stationarity.jl")
for path in ARGS
 started=time();actual=resolve_checkpoint(path);source_sha=open(io->bytes2hex(sha256(io)),actual)
 psi=load_state(path);meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
 family=meta["family"];w=meta["width"];L=meta["length"];ordering=meta["ordering"]
 lat=cylinder(family,L,w;ordering);H=MPO(lat.os,siteinds(psi))
 bonds=unique([max(1,lat.cut-1),lat.cut,min(length(psi)-1,lat.cut+1)])
 records=[finite_bond_stationarity(H,psi,b) for b in bonds]
 @assert open(io->bytes2hex(sha256(io)),actual)==source_sha
 r=Dict("family"=>family,"width"=>w,"length"=>L,"source_checkpoint"=>actual,"source_payload_sha256"=>source_sha,
 "physical_spins"=>length(psi),"bond_dimension"=>maxlinkdim(psi),"number_sector"=>get(meta,"nup",nothing),
 "central_projected_stationarity"=>records,"runtime_seconds"=>time()-started,
 "audit"=>run_provenance(solver="Library ProjMPO saved-state two-site residual",initialization=actual,settings=Dict("sampled_bonds"=>bonds),conserved_quantum_numbers=["Preserved physical U1 sector"]),
 "interpretation"=>"New observable from the same stored state, no optimization. Central retained-basis residuals include truncation and cannot certify global eigenstate, global minimum, or entropy convergence.")
 atomic_json(replace(path,".jls"=>"_stationarity.json"),r)
 println("Saved central stationarity ",path);flush(stdout)
end
