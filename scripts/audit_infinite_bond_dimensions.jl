include("infinite.jl")
records=[]
for file in sort(filter(f->occursin(r"^infinite_.*_chi.*\.json$",basename(f)),readdir("results";join=true)))
 r=JSON3.read(read(file,String),Dict{String,Any});haskey(r,"spatial_entropy") || continue
 get(r,"validation_fixture",false) && continue
 candidate=get(r,"checkpoint_file",replace(file,".json"=>".jls"))
 try
  actual=resolve_checkpoint(candidate);fingerprint=open(io->bytes2hex(sha256(io)),actual)
  psi=load_state(candidate);details=mps_bond_dimension_details(psi)
  @assert open(io->bytes2hex(sha256(io)),actual)==fingerprint
  push!(records,Dict("source_result"=>file,"source_payload"=>actual,"source_payload_sha256"=>fingerprint,
   "family"=>r["family"],"width"=>r["width"],"cap"=>r["cap"],"historical_reported_chi"=>r["chi"],
   "maximum_bond_dimension"=>details["maximum_bond_dimension"],"bond_dimension_details"=>details,"status"=>"audited from unchanged stored tensors"))
 catch exception
  push!(records,Dict("source_result"=>file,"source_payload_requested"=>candidate,"status"=>"not independently recovered","error"=>sprint(showerror,exception)))
 end
end
r=Dict("numerical_optimizations_repeated"=>false,"records"=>records,
 "interpretation"=>"The finite-library helper can omit the periodic wrap bond. Report allocated dimensions for all AL and AR links, separately from physical Schmidt rank and requested caps. Historical raw measurements are retained.",
 "audit"=>run_provenance(solver="Saved-state periodic bond-dimension audit",initialization="Existing infinite checkpoints",settings=Dict(),conserved_quantum_numbers=[]))
atomic_json("results/infinite_bond_dimension_audit.json",r)
println("Audited saved-state bond dimensions: ",count(x->haskey(x,"maximum_bond_dimension"),records),"/",length(records));flush(stdout)
