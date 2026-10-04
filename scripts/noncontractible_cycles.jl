isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/cgs.jl")
# New observable measurements on saved states; no energy optimization is repeated.
for path in ARGS
    started=time();meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    family=meta["family"];L=meta["length"];w=meta["width"];ordering=get(meta,"ordering","axial")
    lat=cylinder(family,L,w;ordering);psi=load_state(path);sites=siteinds(psi)
    settings=get(meta,"solver_settings",get(get(meta,"audit",Dict()),"solver_settings",Dict()))
    basis=get(meta,"basis",get(settings,"basis","physical_spin"))
    if basis=="exact_matter_charge_basis"
        isdefined(Main,:physical_from_charge_basis) || include("../src/matter_charge_basis.jl")
        psi=physical_from_charge_basis(psi,lat);sites=siteinds(psi)
    else
        @assert basis=="physical_spin" "Unrecognized finite-state basis"
    end
    table="geometry/cgs_cycles/$(family)_L$(L)_w$(w)_$(ordering).json"
    definitions=JSON3.read(read(table,String),Dict{String,Any});records=[]
    @assert definitions["family"]==family && definitions["L"]==L && definitions["width"]==w
    for cycle in definitions["cycles"]
        p=cycle["edge_exponents"];gates=cgs_cycle_gates(lat,sites,p)
        rotated=apply(gates,psi;cutoff=1e-13,maxdim=4maxlinkdim(psi))
        z=inner(psi,rotated)/inner(psi,psi)
        probabilities=[real((1+2real(conj(cis(2pi*q/3))*z))/3) for q=0:2]
        norm_error=abs(norm(rotated)-norm(psi))
        push!(records,Dict("label"=>cycle["label"],"winding_number"=>cycle["winding_number"],
            "edge_exponents"=>p,"expectation_real"=>real(z),"expectation_imag"=>imag(z),
            "charge_probabilities"=>probabilities,"rotated_norm_error"=>norm_error,
            "probability_diagnostic_valid"=>norm_error<1e-8 && minimum(probabilities)>-1e-8,
            "charge_purity_error"=>minimum(abs(z-cis(2pi*q/3)) for q=0:2),
            "maxlinkdim_after_apply"=>maxlinkdim(rotated)))
    end
    density=real.(expect(psi,"Sz"))
    nup=hasqns(psi) ? (lat.n+val(flux(psi),"Sz"))÷2 : Int(meta["nup"])
    @assert abs(sum(density)+lat.n/2-nup)<1e-9
    r=Dict("source_result"=>replace(path,".jls"=>".json"),"source_checkpoint"=>resolve_checkpoint(path),
        "family"=>family,"length"=>L,"width"=>w,"ordering"=>ordering,"physical_spins"=>lat.n,
        "actual_nup"=>nup,"bond_dimension"=>maxlinkdim(psi),"geometry_cycle_definitions"=>table,
        "noncontractible_cycle_measurements"=>records,"site_sz_profile_mps_order"=>density,
        "mps_order_physical_spin_ids"=>lat.order,"runtime_seconds"=>time()-started,
        "audit"=>run_provenance(solver="ITensorMPS microscopic circumference-cycle and density contractions",
            settings=Dict("operator_apply_cutoff"=>1e-13,"operator_apply_maxdim"=>4maxlinkdim(psi),"norm_error_target"=>1e-8),
            initialization=path,conserved_quantum_numbers=["total N_up=$nup"]),
        "interpretation"=>"Exact CGS circumference charge diagnostics on a saved variational state; logical/MES identification and global ground selection are not certified.")
    atomic_json(replace(path,".jls"=>"_circumference_cycles.json"),r)
    println(family," L=",L," w=",w," Nup=",nup," circumference cycles=",length(records));flush(stdout)
end
