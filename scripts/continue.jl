isdefined(Main,:cylinder) || include("../src/model.jl")
using Serialization
function continue_cylinder(path,chi)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    family=meta["family"];L=meta["length"];w=meta["width"];seed=get(meta,"seed",7103)
    started=time()
    audit=run_provenance(;seed,solver="ITensorMPS finite DMRG checkpoint continuation",
        settings=Dict("passes"=>2,"sweeps_per_pass"=>8,"cutoff"=>1e-11,"krylovdim"=>12),
        initialization=path,conserved_quantum_numbers=["total Sz"])
    ordering=get(meta,"ordering","axial")
    lat=cylinder(family,L,w;ordering=get(meta,"ordering","axial"));psi=load_state(path)
    H=MPO(lat.os,siteinds(psi))
    records=meta["records"]
    estimate_memory(psi,chi;label="$(family) L$(L) w$(w) continuation")
    for pass=1:2
        before=entropy_at(psi,lat.cut)[1];obs=finite_checkpoint_observer(lat;family,L,w,cap=chi,seed,stage="continuation_pass$(pass)",ordering,audit,cutoff=1e-11,noise=pass==1 ? [1e-7,1e-8,0] : 0)
        e,psi=dmrg(H,psi;nsweeps=8,maxdim=chi,cutoff=1e-11,
            noise=pass==1 ? [1e-7,1e-8,0] : 0,eigsolve_krylovdim=12,outputlevel=1,observer=obs)
        S,p=entropy_at(psi,lat.cut)
        push!(records,Dict("cap"=>chi,"energy"=>e,"entropy"=>S,"maxlinkdim"=>maxlinkdim(psi),
            "schmidt_probabilities"=>p,"sweep_energies"=>energies(obs),"sweep_max_truncation_errors"=>truncerrors(obs),
            "entropy_change_refinement"=>S-before,"process_peak_rss_bytes"=>Sys.maxrss()))
        ordering=get(meta,"ordering","axial")
        output="results/$(family)_L$(L)_w$(w)_chi$(chi)"*(ordering=="axial" ? "" : "_"*ordering)*(get(meta,"seed",7103)==7103 ? "" : "_seed$(meta["seed"])")
        meta["continuation_audit"]=audit;meta["git_commit"]=audit["git_commit"]
        meta["continuation_runtime_seconds"]=time()-started
        open(output*".json","w") do io;JSON3.write(io,meta);end
        manifest=completed_finite_checkpoint(psi,meta,chi,"continuation_pass$(pass)")
        atomic_json(output*".json",manifest)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    continue_cylinder(ARGS[1],parse(Int,ARGS[2]))
end
