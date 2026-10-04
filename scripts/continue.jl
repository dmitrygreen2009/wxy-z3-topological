isdefined(Main,:cylinder) || include("../src/model.jl")
using Serialization
function continue_cylinder(path,chi;cutoff=1e-11,noise_first=[1e-7,1e-8,0],eigsolve_krylovdim=12,eigsolve_maxiter=1,eigsolve_tol=1e-14,eigsolve_verbosity=0)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    family=meta["family"];L=meta["length"];w=meta["width"];seed=get(meta,"seed",7103)
    Random.seed!(seed)
    started=time()
    nup=get(meta,"nup",cld(get(meta,"spins",get(meta,"physical_spins",0)),2))
    audit=run_provenance(;seed,solver="ITensorMPS finite DMRG checkpoint continuation",
        settings=Dict("passes"=>2,"sweeps_per_pass"=>8,"cutoff"=>cutoff,"noise_first_pass"=>noise_first,"krylovdim"=>eigsolve_krylovdim,"eigsolve_maxiter"=>eigsolve_maxiter,"eigsolve_tol"=>eigsolve_tol,"eigsolve_verbosity"=>eigsolve_verbosity,"local_achieved_residual_available"=>false),
        initialization=path,conserved_quantum_numbers=["total N_up=$nup"])
    ordering=get(meta,"ordering","axial")
    lat=cylinder(family,L,w;ordering=get(meta,"ordering","axial"));psi=load_state(path)
    H=MPO(lat.os,siteinds(psi))
    meta["spins"]=lat.n;meta["physical_spins"]=lat.n;meta["vertices"]=lat.nv
    meta["cut"]=lat.cut;meta["legs"]=lat.legs;meta["vertices_table"]=lat.verts;meta["mps_order"]=lat.order
    meta["explicit_number_sector"]=get(meta,"explicit_number_sector",nup!=cld(lat.n,2))
    records=haskey(meta,"records") ? meta["records"] : [Dict("cap"=>meta["cap"],"energy"=>meta["energy"],"entropy"=>meta["entropy"],"maxlinkdim"=>get(meta,"bond_dimension",maxlinkdim(psi)),"schmidt_probabilities"=>get(meta,"schmidt_probabilities",[]),"sweep_energies"=>get(meta,"sweep_energies",[]),"sweep_max_truncation_errors"=>get(meta,"truncation_errors",[]))]
    meta["records"]=records
    estimate_memory(psi,chi;label="$(family) L$(L) w$(w) continuation")
    for pass=1:2
        before=entropy_at(psi,lat.cut)[1];obs=finite_checkpoint_observer(lat;family,L,w,cap=chi,seed,stage="continuation_pass$(pass)",ordering,audit,nup,cutoff,noise=pass==1 ? noise_first : 0)
        e,psi=dmrg(H,psi;nsweeps=8,maxdim=chi,cutoff,
            noise=pass==1 ? noise_first : 0,eigsolve_krylovdim,eigsolve_maxiter,eigsolve_tol,eigsolve_verbosity,outputlevel=1,observer=obs)
        S,p=entropy_at(psi,lat.cut)
        push!(records,Dict("cap"=>chi,"energy"=>e,"entropy"=>S,"maxlinkdim"=>maxlinkdim(psi),
            "schmidt_probabilities"=>p,"sweep_energies"=>energies(obs),"sweep_max_truncation_errors"=>truncerrors(obs),
            "entropy_change_refinement"=>S-before,"process_peak_rss_bytes"=>Sys.maxrss()))
        ordering=get(meta,"ordering","axial")
        output="results/$(family)_L$(L)_w$(w)_chi$(chi)"*(ordering=="axial" ? "" : "_"*ordering)*(get(meta,"seed",7103)==7103 ? "" : "_seed$(meta["seed"])")
        get(meta,"explicit_number_sector",false) && (output*="_Nup$(nup)")
        meta["continuation_audit"]=audit;meta["git_commit"]=audit["git_commit"]
        meta["continuation_runtime_seconds"]=time()-started
        atomic_json(output*".json",meta)
        manifest=completed_finite_checkpoint(psi,meta,chi,"continuation_pass$(pass)")
        atomic_json(output*".json",manifest)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    continue_cylinder(ARGS[1],parse(Int,ARGS[2]))
end
