isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:particle_hole_mps) || include("../src/particle_hole.jl")
function audit_particle_hole_scan(path)
    started=time();scan=JSON3.read(read(path,String),Dict{String,Any})
    reference=JSON3.read(read(replace(scan["reference_file"],".jls"=>".json"),String),Dict{String,Any})
    family=reference["family"];L=reference["length"];w=reference["width"];ordering=get(reference,"ordering","axial")
    lat=cylinder(family,L,w;ordering);records=[]
    for source in scan["records"]
        psi=load_state(source["checkpoint_file"]);mirror=particle_hole_mps(psi)
        H=MPO(lat.os,siteinds(psi));E=real(inner(psi',H,psi));Em=real(inner(mirror',H,mirror))
        S,p=entropy_at(psi,lat.cut);Sm,pm=entropy_at(mirror,lat.cut)
        target=lat.n-source["nup"]
        @assert abs(E-Em)<1e-10 && abs(S-Sm)<1e-10 && abs(norm(mirror)-norm(psi))<1e-10
        @assert val(flux(mirror),"Sz")==2target-lat.n
        audit=run_provenance(;seed=nothing,solver="Exact antiunitary microscopic particle-hole state map",
            settings=Dict("iterations"=>0,"validation_tolerance"=>1e-10),initialization=source["checkpoint_file"],conserved_quantum_numbers=["total N_up=$target"])
        result=Dict("family"=>family,"length"=>L,"width"=>w,"physical_spins"=>lat.n,"ordering"=>ordering,
            "physical_circumference"=>lat.circumference,"nup"=>target,"energy"=>Em,"entropy"=>Sm,"schmidt_probabilities"=>pm,
            "cap"=>source["cap"],"bond_dimension"=>maxlinkdim(mirror),"seed"=>source["seed"],"audit"=>audit,
            "source_checkpoint"=>source["checkpoint_file"],"source_nup"=>source["nup"],"source_audit"=>source["audit"],
            "energy_symmetry_error"=>Em-E,"entropy_symmetry_error"=>Sm-S,"kind"=>"finite")
        checkpoint="results/checkpoints/$(family)_w$(w)_L$(L)_N$(lat.n)_chi$(source["cap"])_Nup$(target)_seed$(source["seed"])_v$(RUN_FORMAT_VERSION)_$(ordering)_particle_hole_partner_complete.jls"
        push!(records,save_checkpoint(checkpoint,mirror,result))
    end
    output=Dict("source_scan"=>path,"records"=>records,"runtime_seconds"=>time()-started,
        "interpretation"=>"Exactly symmetry-related variational states with identical energies and entropies. Independently optimized partner discrepancies reveal solver/sector convergence errors; this map does not certify a global minimum.")
    atomic_json(replace(path,".json"=>"_particle_hole_audit.json"),output)
    println("Verified ",length(records)," exact particle-hole partners for ",path);flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__;for path in ARGS;audit_particle_hole_scan(path);end;end
