include("infinite.jl")
include("../src/canonical_analysis.jl")
for path in ARGS
    started=time();meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    actual=resolve_checkpoint(path);source_sha=open(io->bytes2hex(sha256(io)),actual)
    psi=load_state(path);family=meta["family"];w=meta["width"];cap=meta["cap"]
    tag=get(meta,"measurement_tag","")*"_fixedpoint";ordering=get(meta,"infinite_ordering","matter_first")
    Random.seed!(7114)
    A=hasqns(siteind(psi.AL,1)) ? InfiniteMPS([dense(psi.AL[j]) for j=1:18w],translator(psi.AL)) : psi.AL
    spectrum=transfer_spectrum(A)
    if spectrum["fixed_point_rank_detected"]!=1
        atomic_json("results/infinite_$(family)_w$(w)_chi$(cap)$(tag)_ambiguity.json",
            Dict("source"=>path,"status"=>"nonunique fixed point: preserve stored boundary/sector weights",
                "fixed_point_rank_detected"=>spectrum["fixed_point_rank_detected"],"measurement_git_commit"=>LAUNCH_REVISION))
        continue
    end
    canonical,lambda=canonicalize_left(psi.AL)
    audit=run_provenance(;seed=7114,solver="Official ITensorInfiniteMPS transfer-fixed-point canonicalization",
        settings=Dict("canonicalization_tolerance"=>1e-14),initialization=path,conserved_quantum_numbers=hasqns(siteind(psi.AL,1)) ? ["total Sz"] : String[])
    audit["number_background"]=get(meta,"number_background",infinite_number_background(psi))
    audit["source_payload"]=actual;audit["source_payload_sha256"]=source_sha
    audit["source_metadata_sha256"]=bytes2hex(sha256(read(replace(path,".jls"=>".json"))))
    audit["optimization_u1_conserving_ansatz"]=hasqns(siteind(psi.AL,1))
    audit["measurement_densification_only"]=hasqns(siteind(psi.AL,1))
    audit["source_variational_audit"]=get(meta,"audit",nothing)
    audit["source_spatial_entropy"]=meta["spatial_entropy"]
    audit["source_energy_cell"]=meta["energy_cell"]
    audit["source_canonical_error"]=meta["canonical_error"]
    audit["canonical_solver_residual"]=get(meta,"solver_residual",get(meta,"canonical_error",nothing))
    audit["runtime_seconds"]=time()-started
    audit["interpretation"]="No variational improvement. Recanonicalize the same AL state; inherited VUMPS residual still controls variational convergence. Nonunique fixed points are skipped to avoid changing sector weights."
    atomic_json("results/infinite_$(family)_w$(w)_chi$(cap)$(tag)_solver_audit.json",audit)
    H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering),siteinds(only,canonical))
    @assert open(io->bytes2hex(sha256(io)),actual)==source_sha "Source checkpoint changed during canonicalization"
    measure_infinite(canonical,H,family,w,cap,meta["iteration"];tag,ordering)
end
