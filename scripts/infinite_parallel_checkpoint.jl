include("infinite.jl")
for path in ARGS
    actual=resolve_checkpoint(path);fingerprint=open(io->bytes2hex(sha256(io)),actual)
    psi=load_state(path);meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    family=meta["family"];w=meta["width"];cap=meta["cap"];ordering=get(meta,"infinite_ordering","star")
    tag=replace(get(meta,"measurement_tag",""),"_fixedpoint"=>"")*"_parallel"
    seed=7232;Random.seed!(seed)
    estimate_memory(psi,cap;label="Parallel VUMPS continuation from immutable checkpoint")
    H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering,cell_slices=length(psi.AL)÷(9w)),siteinds(only,psi))
    audit=run_provenance(seed=seed,solver="Official ITensorInfiniteMPS parallel VUMPS continuation",
        settings=Dict("tol"=>1e-7,"maxiter"=>40,"local_eigensolver_tolerance"=>1e-10,"multisite_update_algorithm"=>"parallel","bond_dimension_expanded"=>false),
        initialization=actual,conserved_quantum_numbers=hasqns(siteind(psi.AL,1)) ? ["U1 mean physical density"] : String[])
    audit["source_payload_sha256"]=fingerprint;audit["source_solver_residual"]=get(meta,"solver_residual",nothing)
    audit["interpretation"]="Independent algorithm branch from the existing checkpoint; source results preserved. Simultaneous tensor updates share library environments, unlike sequential updates. No tolerance is loosened."
    psi=audited_vumps(H,psi;family,w,cap,maxiter=40,tol=1e-7,solver_tol=x->1e-10,ordering,tag,seed,audit,update_algorithm="parallel")
    @assert open(io->bytes2hex(sha256(io)),actual)==fingerprint
    measure_infinite(psi,H,family,w,cap,1;ordering,tag)
end
