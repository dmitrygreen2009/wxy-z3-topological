# Convert and project a saved physical MPS; never modify its source payload.
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/winding_checkpoint_bridge.jl")
function convert_winding_checkpoint(path,charge;maxdim=nothing)
    Random.seed!(7256);started=time();actual=resolve_checkpoint(path);source_sha=open(io->bytes2hex(sha256(io)),actual)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    settings=get(meta,"solver_settings",get(get(meta,"audit",Dict()),"solver_settings",Dict()))
    @assert get(meta,"basis",get(settings,"basis","physical_spin"))=="physical_spin"
    family=meta["family"];L=meta["length"];w=meta["width"];ordering=get(meta,"ordering","axial")
    @assert ordering=="star" "Direct sector continuation currently uses audited star ordering"
    lat=cylinder(family,L,w;ordering);psi=load_state(path);nup=(lat.n+val(flux(psi),"Sz"))÷2
    table="geometry/cgs_cycles/$(family)_L$(L)_w$(w)_$(ordering).json"
    defs=JSON3.read(read(table,String),Dict{String,Any});cycle=first(c for c in defs["cycles"] if c["winding_number"]==1)
    p=Int.(cycle["edge_exponents"]);cap=maxdim===nothing ? max(64,4maxlinkdim(psi)) : maxdim
    estimate_memory(psi,cap;label="Checkpoint number/winding projection and exact matter basis bridge")
    Ebefore=real(inner(psi',MPO(lat.os,siteinds(psi)),psi));Sbefore,_=entropy_at(deepcopy(psi),lat.cut)
    converted,record=bridge_physical_checkpoint(psi,lat,p;charge,maxdim=cap,cutoff=1e-13)
    H=MPO(charge_basis_opsum(lat),siteinds(converted));Eafter=real(inner(converted',H,converted))
    S,probs=entropy_at(converted,lat.cut)
    if abs(record["projection_weight"]-1)<1e-10
        @assert abs(Eafter-Ebefore)<1e-9 && abs(S-Sbefore)<1e-9 "A pure-sector basis conversion changed physical energy/entropy"
    end
    @assert open(io->bytes2hex(sha256(io)),actual)==source_sha
    seed=7256
    audit=run_provenance(;seed,solver="ITensorMPS exact P_q B† checkpoint conversion",
        settings=Dict("basis"=>"exact_matter_charge_basis","winding_charge"=>charge,"onsite_winding_weights"=>record["onsite_winding_weights"],
            "edge_exponents"=>p,"cutoff"=>1e-13,"conversion_maxdim"=>cap,"krylovdim"=>12,"eigsolve_maxiter"=>30,"eigsolve_tol"=>1e-13),
        initialization=actual,conserved_quantum_numbers=["Physical N_up=$nup","Microscopic winding charge=$charge"])
    audit["source_payload_sha256"]=source_sha
    result=Dict("family"=>family,"length"=>L,"width"=>w,"physical_circumference"=>lat.circumference,"spins"=>lat.n,"vertices"=>lat.nv,
        "ordering"=>ordering,"seed"=>seed,"basis"=>"exact_matter_charge_basis","nup"=>nup,"total_sz"=>nup-lat.n/2,"filling_fraction"=>nup/lat.n,
        "source_payload"=>actual,"source_payload_sha256"=>source_sha,"source_energy"=>Ebefore,"source_entropy"=>Sbefore,
        "conversion_record"=>record,"conversion_energy"=>Eafter,"conversion_entropy"=>S,"schmidt_probabilities"=>probs,
        "geometry_cycle_definitions"=>table,"loop_label"=>cycle["label"],"audit"=>audit,"git_commit"=>audit["git_commit"],
        "runtime_seconds"=>time()-started,"interpretation"=>"New sector-projected branch from an unchanged original checkpoint. No optimization or convergence certificate.")
    manifest=completed_finite_checkpoint(converted,result,cap,"winding$(charge)_basis_bridge")
    output="results/$(family)_L$(L)_w$(w)_Nup$(nup)_winding$(charge)_basis_bridge.json"
    atomic_json(output,manifest);println("Converted checkpoint: ",manifest["checkpoint_file"]);flush(stdout)
    manifest["checkpoint_file"]
end
if abspath(PROGRAM_FILE)==@__FILE__
    @assert length(ARGS)==2
    convert_winding_checkpoint(ARGS[1],parse(Int,ARGS[2]))
end
