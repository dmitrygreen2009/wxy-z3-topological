include("infinite.jl")
include("../src/infinite_observables.jl")
# New two-point contractions, using previously measured one-point data.
path=ARGS[1];onepoint_path=ARGS[2];started=time()
meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
@assert get(meta,"infinite_ordering","matter_first")=="star"
old=JSON3.read(read(onepoint_path,String),Dict{String,Any})
if haskey(old,"source_payload_sha256")
    @assert open(io->bytes2hex(sha256(io)),resolve_checkpoint(path))==old["source_payload_sha256"] "Stored one-point source differs from this checkpoint"
end

@assert old["source_file"]==path "One-point data must describe this same saved state"
psi=load_state(path);psi=InfiniteCanonicalMPS(dense(psi.AL),dense(psi.C),dense(psi.AR))
sites=siteinds(only,psi);w=meta["width"];records=[]
for (k,first) in enumerate([1,6])
    onepoint=complex(old["matter_triplet_expectations_real"][k],old["matter_triplet_expectations_imag"][k])
    for cells in [1,2,4,8]
        raw=expect(psi,infinite_triplet_pair_mpo(sites,first,first+length(psi.AL)*cells))
        connected=raw-abs2(onepoint)
        push!(records,Dict("first_matter_site"=>first,"distance_cells"=>cells,
            "one_point_real"=>real(onepoint),"one_point_imag"=>imag(onepoint),
            "correlation_real"=>real(raw),"correlation_imag"=>imag(raw),
            "connected_real"=>real(connected),"connected_imag"=>imag(connected)))
    end
end
actual=resolve_checkpoint(path)
result=Dict("family"=>meta["family"],"width"=>w,"cell_slices"=>length(psi.AL)÷(9w),"source_checkpoint"=>actual,
    "source_payload_sha256"=>open(io->bytes2hex(sha256(io)),actual),
    "previous_one_point_data"=>onepoint_path,
    "previous_one_point_payload_fingerprint_available"=>haskey(old,"source_payload_sha256"),"previous_one_point_sha256"=>bytes2hex(sha256(read(onepoint_path))),
    "source_canonical_error"=>meta["canonical_error"],"source_solver_residual"=>get(meta,"solver_residual",nothing),
    "matter_charge_three_correlations"=>records,"runtime_seconds"=>time()-started,
    "audit"=>run_provenance(solver="ITensorInfiniteMPS CGS-invariant charge-three two-point contractions",
        settings=Dict("distances_cells"=>[1,2,4,8]),initialization=path,conserved_quantum_numbers=[]),
    "interpretation"=>"Physical invariant connected correlations of this variational state; compared with its full transfer length without interpreting either as a 2D spectral gap. Previously validated one-point contractions are reused.")
atomic_json(replace(path,".jls"=>"_physical_correlations.json"),result)
println("Saved invariant two-point correlations for ",path);flush(stdout)
