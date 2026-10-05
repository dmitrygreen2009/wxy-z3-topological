# Saved-state measurement only. Never repeat an energy optimization.
include("infinite.jl")
include("../src/cgs_operator_mpo.jl")
include("../src/canonical_analysis.jl")
include("../src/transfer_uniqueness.jl")
function reaudit_infinite_loops(path)
    started=time();meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    @assert checkpoint_basis(meta)=="physical_spin" "Use the winding-basis measurement driver for rotated infinite states"
    @assert get(meta,"infinite_ordering","matter_first")=="star"
    actual=resolve_checkpoint(path);source_sha=open(io->bytes2hex(sha256(io)),actual)
    psi=load_state(path);family=meta["family"];w=meta["width"];n=length(psi.AL);slices=n÷(9w)
    background=get(meta,"number_background",infinite_number_background(psi))
    canonical_status="Stored state retained"
    if meta["canonical_error"]>1e-10
        A=hasqns(siteind(psi.AL,1)) ? InfiniteMPS([dense(psi.AL[j]) for j=1:n],translator(psi.AL)) : psi.AL
        certificate=transfer_uniqueness_certificate(A)
        if !certificate["primitive_peripheral_spectrum_certified"]
            atomic_json(replace(path,".jls"=>"_loop_reaudit.json"),Dict("source"=>path,"source_payload_sha256"=>source_sha,
                "status"=>"Skipped: primitive fixed point not certified; no arbitrary sector/boundary change","transfer_uniqueness_certificate"=>certificate))
            return
        end
        psi,_=canonicalize_left(psi.AL);canonical_status="Same AL state recanonicalized; no optimization"
    else
        psi=InfiniteCanonicalMPS(dense(psi.AL),dense(psi.C),dense(psi.AR))
    end
    sites=siteinds(only,psi);records=[]
    sz_profile=[real(expect(psi,"Sz",j)) for j=1:n]
    total_sz=sum(sz_profile);cell_nup=n/2+total_sz
    # Include one complete translated cell to expose finite-index aliases.
    for offset in 0:9w:n
        spec=cgs_circumference_spec(family,w;offset)
        z=expect(psi,cgs_product_mpo(sites,spec))
        push!(records,Dict("offset_sites"=>offset,"support_first"=>spec.first_site,"support_last"=>spec.last_site,
            "matter_transformations"=>spec.triplets,"gauge_exponents"=>spec.gauges,"winding_number"=>1,
            "expectation_real"=>real(z),"expectation_imag"=>imag(z),"unitary_variance"=>1-abs2(z),
            "charge_probabilities"=>[(1+2real(conj(cis(2pi*q/3))*z))/3 for q=0:2],
            "purity_error"=>minimum(abs(z-cis(2pi*q/3)) for q=0:2)))
    end
    @assert abs(complex(records[1]["expectation_real"],records[1]["expectation_imag"])-complex(records[end]["expectation_real"],records[end]["expectation_imag"]))<1e-10 "Full-cell translation changed loop expectation"
    plaquettes=[]
    table="geometry/cgs_plaquettes/$(family)_L4_w$(w)_star.json"
    if isfile(table)
        definitions=JSON3.read(read(table,String),Dict{String,Any})
        for cycle in definitions["plaquettes"]
            haskey(cycle,"infinite_star_operator") || continue
            r=cycle["infinite_star_operator"]
            spec=(;triplets=[Tuple(Int.(x)) for x in r["triplets"]],gauges=[Tuple(Int.(x)) for x in r["gauges"]],first_site=r["first_site"],last_site=r["last_site"])
            z=expect(psi,cgs_product_mpo(sites,spec))
            push!(plaquettes,Dict("label"=>cycle["label"],"winding_number"=>0,"expectation_real"=>real(z),"expectation_imag"=>imag(z),
                "unitary_variance"=>1-abs2(z),"charge_probabilities"=>[(1+2real(conj(cis(2pi*q/3))*z))/3 for q=0:2]))
        end
    end
    pairs=[]
    for cells in [1,2,4,8,16]
        firstspec=cgs_circumference_spec(family,w)
        secondspec=cgs_circumference_spec(family,w;offset=cells*n)
        # Armchair loops occupy two slices; adjacent one-slice cells overlap.
        secondspec.first_site>firstspec.last_site || continue
        z=expect(psi,cgs_product_mpo(sites,cgs_disjoint_pair_spec(firstspec,secondspec)))
        mean=complex(records[1]["expectation_real"],records[1]["expectation_imag"])
        connected=z-abs2(mean)
        push!(pairs,Dict("distance_cells"=>cells,"distance_slices"=>cells*slices,
            "real"=>real(z),"imag"=>imag(z),"connected_real"=>real(connected),"connected_imag"=>imag(connected)))
    end
    @assert open(io->bytes2hex(sha256(io)),actual)==source_sha "Source checkpoint changed: discard measurement"
    result=Dict("source"=>path,"source_payload"=>actual,"source_payload_sha256"=>source_sha,"family"=>family,"width"=>w,
        "cell_slices"=>slices,"cell_spins"=>n,"canonical_action"=>canonical_status,
        "source_canonical_error"=>meta["canonical_error"],"source_solver_residual"=>get(meta,"solver_residual",nothing),
        "winding_measurements"=>records,"contractible_measurements"=>plaquettes,"winding_pair_correlations"=>pairs,
        "number_background"=>background,"unit_cell_nup_expectation"=>cell_nup,
        "unit_cell_sz_expectation"=>total_sz,"filling_fraction"=>cell_nup/n,
        "filling_selection"=>get(background,"density_enforced",nothing)===true ? "Fixed rational U1 MPS background" : get(meta,"u1_conserving_ansatz",false) ? "Legacy U1 background unresolved" : "Variational number expectation, not a fixed QN sector",
        "site_sz_profile"=>sz_profile,"runtime_seconds"=>time()-started,
        "audit"=>run_provenance(solver="Official InfiniteMPS contractions with scalar translated indices",settings=Dict("operator_cutoff"=>1e-14),
            initialization=path,conserved_quantum_numbers=[]),
        "interpretation"=>"Remeasured exact microscopic CGS operators on the unchanged saved variational state. Sector/MES and global filling identification remain separate gates. Archived triplets excluded.")
    atomic_json(replace(path,".jls"=>"_loop_reaudit.json"),result)
    println("Loop re-audit saved: ",path);flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__
    for path in ARGS;reaudit_infinite_loops(path);GC.gc();end
end
