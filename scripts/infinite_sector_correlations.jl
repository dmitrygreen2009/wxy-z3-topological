include("infinite.jl")
include("../src/cgs_operator_mpo.jl")
# Diagnose whether the narrow-quotient transfer length couples to loop sectors.
path=ARGS[1];onepoint_path=ARGS[2];started=time()
meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
@assert meta["family"]=="zigzag" && meta["width"]==1 && meta["infinite_ordering"]=="star"
old=JSON3.read(read(onepoint_path,String),Dict{String,Any})
if haskey(old,"source_payload_sha256")
    @assert open(io->bytes2hex(sha256(io)),resolve_checkpoint(path))==old["source_payload_sha256"] "Stored one-point source differs from this checkpoint"
end
@assert old["source_file"]==path
z=complex(old["parallel_two_cycle_expectations_real"][1],old["parallel_two_cycle_expectations_imag"][1])
actual=resolve_checkpoint(path);psi=load_state(path)
psi=InfiniteCanonicalMPS(dense(psi.AL),dense(psi.C),dense(psi.AR));sites=siteinds(only,psi)
records=[]
for cells in [1,2,4,8,16]
    spec=cgs_disjoint_pair_spec(cgs_circumference_spec("zigzag",1),cgs_circumference_spec("zigzag",1;offset=18cells))
    raw=expect(psi,cgs_product_mpo(sites,spec));connected=raw-abs2(z)
    push!(records,Dict("distance_cells"=>cells,"distance_slices"=>2cells,
        "correlation_real"=>real(raw),"correlation_imag"=>imag(raw),
        "connected_real"=>real(connected),"connected_imag"=>imag(connected)))
end
result=Dict("family"=>"zigzag","width"=>1,"source_checkpoint"=>actual,
    "source_payload_sha256"=>open(io->bytes2hex(sha256(io)),actual),"previous_one_point_data"=>onepoint_path,
    "previous_one_point_payload_fingerprint_available"=>haskey(old,"source_payload_sha256"),
    "previous_one_point_json_sha256"=>bytes2hex(sha256(read(onepoint_path))),
    "source_canonical_error"=>meta["canonical_error"],"source_solver_residual"=>get(meta,"solver_residual",nothing),
    "source_xi_cells"=>meta["xi_cells"],"loop_connected_correlations"=>records,
    "runtime_seconds"=>time()-started,"audit"=>run_provenance(solver="ITensorInfiniteMPS winding-loop adjoint-pair contractions",
        settings=Dict("operator_cutoff"=>1e-14,"distances_cells"=>[1,2,4,8,16]),initialization=path,conserved_quantum_numbers=[]),
    "interpretation"=>"New connected loop-sector correlations on the saved narrow-quotient variational state, reusing its existing one-point contraction. Winding operators have local support only on this special width-one quotient. No 2D gap, order, or topology conclusion follows.")
atomic_json(replace(path,".jls"=>"_sector_correlations.json"),result)
println("Saved loop-sector correlations for ",path);flush(stdout)
