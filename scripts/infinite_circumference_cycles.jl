include("infinite.jl")
include("../src/cgs_operator_mpo.jl")
include("../src/infinite_observables.jl")
# Measurement only: preserve optimized checkpoints and their solver residuals.
for path in ARGS
    started=time()
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    @assert get(meta,"infinite_ordering","matter_first")=="star" "Operator coordinates require the explicitly audited star ordering"
    actual=resolve_checkpoint(path);psi=load_state(path)
    psi=InfiniteCanonicalMPS(dense(psi.AL),dense(psi.C),dense(psi.AR))
    sites=siteinds(only,psi);w=meta["width"];family=meta["family"]
    records=[]
    for offset in [0,9w]
        spec=cgs_circumference_spec(family,w;offset)
        operator=cgs_product_mpo(sites,spec)
        z=expect(psi,operator)
        probabilities=[(1+2real(conj(cis(2pi*q/3))*z))/3 for q=0:2]
        push!(records,Dict("offset_sites"=>offset,"winding_number"=>1,
            "matter_transformations"=>spec.triplets,"gauge_exponents"=>spec.gauges,
            "expectation_real"=>real(z),"expectation_imag"=>imag(z),
            "charge_probabilities"=>probabilities,"operator_max_bond_dimension"=>maxlinkdim(operator),
            "charge_purity_error"=>minimum(abs(z-cis(2pi*q/3)) for q=0:2)))
    end
    result=Dict("family"=>family,"width"=>w,"source_checkpoint"=>path,
        "source_payload"=>actual,"source_sha256"=>open(io->bytes2hex(SHA.sha256(io)),actual),
        "source_canonical_error"=>meta["canonical_error"],
        "source_solver_residual"=>get(meta,"solver_residual",nothing),
        "noncontractible_cycle_measurements"=>records,
        "runtime_seconds"=>time()-started,
        "audit"=>run_provenance(solver="ITensorInfiniteMPS factorized microscopic circumference MPO contraction",
            settings=Dict("operator_cutoff"=>1e-14),initialization=path,conserved_quantum_numbers=[]),
        "interpretation"=>"Circumference CGS charges of this saved variational state. Values remain provisional when canonical or solver residuals are large. A loop eigenstate alone does not certify a minimally entangled state or a 2D phase.")
    atomic_json(replace(path,".jls"=>"_circumference_cycles.json"),result)
    println(family," w=",w," circumference charge diagnostics saved");flush(stdout)
end
