# Energy remeasurement required by a failed prior canonical-consistency gate.
# No optimizer; same physical AL amplitudes and certified primitive fixed point.
include("infinite.jl")
include("../src/canonical_analysis.jl")
include("../src/transfer_uniqueness.jl")
function energy48()
 started=time();path="results/infinite_armchair_w2_chi64_qn_star_active.jls";meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
 @assert !meta["bound_comparison_state_consistency_checked"]
 actual=resolve_checkpoint(path);digest=open(io->bytes2hex(sha256(io)),actual)
 prior=JSON3.read(read("results/infinite_armchair_w2_chi64_qn_star_active_loop_reaudit.json",String),Dict{String,Any})
 @assert prior["source_payload_sha256"]==digest && prior["transfer_uniqueness_certificate"]["primitive_peripheral_spectrum_certified"]
 source=load_state(path);psi,_=canonicalize_left(source.AL);n=length(psi.AL);slices=n÷(9meta["width"])
 H=InfiniteSum{MPO}(infinite_opsum(meta["family"],meta["width"];ordering=meta["infinite_ordering"],cell_slices=slices),siteinds(only,psi))
 energies=expect(psi,H);energy=real(sum(energies));canon=Float64[];left=Float64[];right=Float64[];centers=Float64[]
 for j=1:n
  A=psi.AL[j]*psi.C[j];B=psi.C[j-1]*psi.AR[j];z=inner(A,B);push!(canon,norm(A-(abs(z)>0 ? conj(z)/abs(z) : 1)*B))
  r=commonind(psi.AL[j],psi.AL[j+1]);l=commonind(psi.AR[j-1],psi.AR[j])
  push!(left,norm(dense(psi.AL[j]*dag(prime(psi.AL[j],r)))-dense(delta(r,dag(prime(r)))))/sqrt(dim(r)))
  push!(right,norm(dense(psi.AR[j]*dag(prime(psi.AR[j],l)))-dense(delta(l,dag(prime(l)))))/sqrt(dim(l)))
  push!(centers,abs(norm(psi.C[j])-1))
 end
 cert=transfer_uniqueness_certificate(psi.AL);tnorm=abs(hypot(cert["values_real"][1],cert["values_imag"][1])-1)
 errors=Dict("canonical_error"=>maximum(canon),"left_isometry_error"=>maximum(left),"right_isometry_error"=>maximum(right),"center_normalization_error"=>maximum(centers),"transfer_normalization_error"=>tnorm)
 @assert cert["primitive_peripheral_spectrum_certified"] && maximum(values(errors))<1e-10
 sz=sum(real(expect(psi,"Sz",j)) for j=1:n);@assert abs(sz)<1e-10
 pervertex=energy/(2slices*meta["width"]);trial=-1.5
 @assert open(io->bytes2hex(sha256(io)),actual)==digest
 result=Dict("source"=>path,"source_payload_sha256"=>digest,"reason"=>"Prior bound-comparison gate was false at canonical error 3.4e-9; new exact coherent trial warrants consistent remeasurement, not an optimization", "canonical_action"=>"Official same-AL recanonicalization", "energy_cell"=>energy,"energy_per_vertex"=>pervertex,"energy_change_vs_raw_saved"=>energy-meta["energy_cell"],"state_consistency_errors"=>errors,"unchanged_consistency_tolerance"=>1e-10,"transfer_certificate"=>cert,"cell_spins"=>n,"unit_cell_nup_expectation"=>n/2+sz,"unit_cell_sz_expectation"=>sz,"filling_fraction"=>.5+sz/n,"filling_selection"=>"Original fixed rational U1 background; same physical state", "same_cylinder_coherent_trial_energy_per_vertex"=>trial,"energy_above_trial_per_vertex"=>pervertex-trial,"excluded_as_unrestricted_global_ground_candidate"=>pervertex>trial+1e-8,
 "qualification"=>"Comparison is with a valid number-projected coherent trial at half filling on the same infinite cylinder, over all CGS sectors. It does not rule out stationarity or a sector-restricted minimum in the candidate's particular cycle background. No 2D tiled bound imported into the cylinder.","optimization_performed"=>false,"runtime_seconds"=>time()-started,"git_commit"=>LAUNCH_REVISION)
 atomic_json("results/physics48_infinite_energy.json",result);println("Consistent saved candidate energy / vertex ",pervertex," exceeds half-filled trial by ",pervertex-trial)
end
if abspath(PROGRAM_FILE)==@__FILE__;energy48();end
