include("infinite.jl")
function parallel_refine_checkpoint(path,target=64)
actual=resolve_checkpoint(path);source_sha=open(io->bytes2hex(sha256(io)),actual)
psi=load_state(path);meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
family=meta["family"];w=meta["width"];ordering=get(meta,"infinite_ordering","star");slices=length(psi.AL)÷(9w)
tag=get(meta,"measurement_tag","")*"_expanded"
H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering,cell_slices=slices),siteinds(only,psi))
caps=sort(unique(vcat(meta["cap"],[c for c in [16,32,64,128,256] if meta["cap"]<c<=target],target)))
@assert minimum(caps)>=meta["cap"]
for (stage,cap) in enumerate(caps)
 seed=7240+stage;Random.seed!(seed)
 estimate_memory(psi,cap;label="Period-aware official parallel VUMPS refinement")
 psi=subspace_expansion(psi,H;cutoff=1e-10,maxdim=cap)
 audit=run_provenance(seed=seed,solver="Official ITensorInfiniteMPS parallel VUMPS expanded continuation",
  settings=Dict("tol"=>1e-7,"maxiter"=>40,"local_eigensolver_tolerance"=>1e-10,"subspace_expansion_cutoff"=>1e-10,"cap"=>cap,"cell_slices"=>slices),
  initialization=actual,conserved_quantum_numbers=hasqns(siteind(psi.AL,1)) ? ["U1 physical density"] : String[])
 audit["source_payload_sha256"]=source_sha
 audit["interpretation"]="Additional sweeps and bond space from preserved checkpoint; canonical and achieved local tolerances unchanged. Independent branch, no ground-state or MES certificate."
 psi=audited_vumps(H,psi;family,w,cap,ordering,tag,seed,audit,maxiter=40,tol=1e-7,solver_tol=x->1e-10,update_algorithm="parallel")
 @assert open(io->bytes2hex(sha256(io)),actual)==source_sha
 measure_infinite(psi,H,family,w,cap,stage;ordering,tag)
end

end
if abspath(PROGRAM_FILE)==@__FILE__
 parallel_refine_checkpoint(ARGS[1],length(ARGS)>=2 ? parse(Int,ARGS[2]) : 64)
end
