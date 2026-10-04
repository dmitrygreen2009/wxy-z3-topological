include("infinite.jl")
include("../src/cgs_operator_mpo.jl")
# New measurements only; no optimization, entropy, or transfer spectrum rerun.
for path in ARGS
    started=time();meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    @assert get(meta,"infinite_ordering","matter_first")=="star"
    family=meta["family"];w=meta["width"];actual=resolve_checkpoint(path)
    psi=load_state(path);psi=InfiniteCanonicalMPS(dense(psi.AL),dense(psi.C),dense(psi.AR))
    sites=siteinds(only,psi)
    definition="geometry/cgs_plaquettes/$(family)_L4_w$(w)_star.json"
    table=JSON3.read(read(definition,String),Dict{String,Any});records=[]
    for cycle in table["plaquettes"]
        haskey(cycle,"infinite_star_operator") || continue
        d=cycle["infinite_star_operator"]
        spec=(triplets=[Tuple(Int.(x)) for x in d["triplets"]],gauges=[Tuple(Int.(x)) for x in d["gauges"]],first_site=d["first_site"],last_site=d["last_site"])
        operator=cgs_product_mpo(sites,spec);z=expect(psi,operator)
        push!(records,Dict("label"=>cycle["label"],"winding_number"=>0,
            "expectation_real"=>real(z),"expectation_imag"=>imag(z),
            "charge_probabilities"=>[(1+2real(conj(cis(2pi*q/3))*z))/3 for q=0:2],
            "charge_purity_error"=>minimum(abs(z-cis(2pi*q/3)) for q=0:2)))
    end
    result=Dict("family"=>family,"width"=>w,"source_result"=>replace(path,".jls"=>".json"),
        "source_checkpoint"=>actual,"source_payload_sha256"=>open(io->bytes2hex(sha256(io)),actual),
        "geometry_definitions"=>definition,"geometry_sha256"=>bytes2hex(sha256(read(definition))),
        "source_canonical_error"=>meta["canonical_error"],"source_solver_residual"=>get(meta,"solver_residual",nothing),
        "contractible_cycle_measurements"=>records,"runtime_seconds"=>time()-started,
        "audit"=>run_provenance(solver="ITensorInfiniteMPS microscopic contractible CGS operator contractions",
            settings=Dict("operator_cutoff"=>1e-14),initialization=path,conserved_quantum_numbers=[]),
        "interpretation"=>"Contractible CGS charge weights of the saved variational state; no assumption that the lowest physical sector has charge zero, and no thermodynamic or topological phase certificate.")
    atomic_json(replace(path,".jls"=>"_contractible_cycles.json"),result)
    println(family," w=",w," contractible cycles measured=",length(records));flush(stdout)
end
