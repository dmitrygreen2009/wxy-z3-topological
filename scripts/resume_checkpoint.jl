include("cylinders.jl")
include("continue.jl")
include("infinite.jl")
include("../src/resume_controls.jl")
function newest_checkpoint(family,L,w;nup=nothing)
    candidates=[]
    for file in readdir("results/checkpoints";join=true)
        endswith(file,".json") || endswith(file,".json.previous") || continue
        meta=try JSON3.read(read(file,String),Dict{String,Any}) catch; continue end
        get(meta,"family",nothing)==family && get(meta,"length",nothing)==L && get(meta,"width",nothing)==w || continue
        get(meta,"kind",nothing)=="finite" || continue
        neutral=cld(get(meta,"physical_spins",get(meta,"spins",0)),2)
        requested=nup===nothing ? neutral : nup
        get(meta,"nup",neutral)==requested || continue
        path=replace(file,".json"=>".jls");isfile(path) && push!(candidates,path)
    end
    for path in sort(candidates;by=mtime,rev=true)
        try load_valid_checkpoint(path);return path catch error
            println("Skipping invalid checkpoint ",path,": ",sprint(showerror,error))
        end
    end
    error("No valid checkpoint for $family L=$L w=$w")
end
function resume_checkpoint(path,target=nothing)
    psi,meta=load_valid_checkpoint(path)
    family=meta["family"];w=meta["width"];cap=meta["cap"]
    if meta["kind"]=="infinite"
        ordering=meta["infinite_ordering"];tag=get(meta,"measurement_tag","")
        H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering,cell_slices=length(psi.AL)÷(9w)),siteinds(only,psi))
        settings=get(meta,"solver_settings",Dict())
        rule=get(settings,"local_solver_tolerance_rule","historical: use fixed 1e-10")
        algorithm=get(settings,"multisite_update_algorithm","sequential")
        local_tol=min(get(settings,"local_eigensolver_tolerance_at_initial_residual",1e-10),1e-10)
        local_solver=rule=="max(residual/10,1e-10)" ? (x->max(x/10,1e-10)) : (x->local_tol)
        seed=get(meta,"random_seed",get(meta,"seed",7103));Random.seed!(seed)
        remaining=max(1,get(meta,"stage_maxiter",40)-meta["iteration"])
        function resumed_audit(nextcap,iterations)
            run_provenance(;seed,solver="ITensorInfiniteMPS VUMPS checkpoint resume",
                settings=merge(copy(settings),Dict("maxiter"=>iterations,"cap"=>nextcap,"source_iteration"=>meta["iteration"],"local_solver_tolerance_rule"=>rule,"multisite_update_algorithm"=>algorithm)),
                initialization=path,conserved_quantum_numbers=hasqns(siteind(psi.AL,1)) ? ["U1 background retained from checkpoint"] : String[])
        end
        audit=resumed_audit(cap,remaining)
        audit["source_checkpoint_sha256"]=get(meta,"checkpoint_sha256",nothing)
        audit["source_git_commit"]=get(meta,"git_commit",nothing)
        psi=audited_vumps(H,psi;family,w,cap,maxiter=remaining,seed,audit,
            tol=get(meta,"tol",1e-7),solver_tol=local_solver,solver_tolerance_rule=rule,ordering,tag,update_algorithm=algorithm)
        measure_infinite(psi,H,family,w,cap,meta["iteration"]+1;ordering,tag)
        target===nothing && return
        for nextcap in [16,32,64,128,256,512]
            cap<nextcap<=target || continue
            estimate_memory(psi,nextcap;label="$(family) infinite resumed w$(w) expansion")
            psi=subspace_expansion(psi,H;cutoff=min(get(settings,"subspace_expansion_cutoff",1e-10),1e-10),maxdim=nextcap)
            psi=audited_vumps(H,psi;family,w,cap=nextcap,maxiter=40,seed,audit=resumed_audit(nextcap,40),tol=get(meta,"tol",1e-7),solver_tol=local_solver,solver_tolerance_rule=rule,ordering,tag,update_algorithm=algorithm)
            measure_infinite(psi,H,family,w,nextcap,100+nextcap;ordering,tag)
        end
        return
    end
    L=meta["length"];nup=get(meta,"nup",cld(meta["physical_spins"],2));explicit_sector=get(meta,"explicit_number_sector",nup!=cld(meta["physical_spins"],2));seed=get(meta,"random_seed",get(meta,"seed",7103));ordering=meta["ordering"]
    Random.seed!(seed)
    lat=cylinder(family,L,w;ordering);H=MPO(lat.os,siteinds(psi))
    phase=string(get(meta,"stage",get(meta,"phase","resume")))
    settings=get(meta,"solver_settings",Dict())
    planned=occursin("refinement",phase) ? 6 : get(settings,"sweeps_per_cap",get(settings,"stage_sweeps",get(settings,"sweeps_per_pass",get(settings,"sweeps",8))))
    remaining=max(1,planned-get(meta,"sweep",0))
    target===nothing && (target=get(settings,"target_chi",get(settings,"chi",cap)))
    controls=finite_resume_controls(settings,phase)
    audit=run_provenance(;seed,solver="ITensorMPS DMRG resume",settings=merge(copy(settings),Dict("remaining_sweeps"=>remaining,"cutoff"=>get(meta,"cutoff",1e-11),"target_chi"=>target,"krylovdim"=>controls.eigsolve_krylovdim,"eigsolve_maxiter"=>controls.eigsolve_maxiter,"eigsolve_tol"=>controls.eigsolve_tol,"eigsolve_verbosity"=>controls.eigsolve_verbosity)),
        initialization=path,conserved_quantum_numbers=["total N_up=$nup"])
    cutoff=get(meta,"cutoff",1e-11)
    obs=finite_checkpoint_observer(lat;family,L,w,cap,seed,stage="resumed",ordering,audit,nup,cutoff,noise=0)
    started=time()
    E,psi=dmrg(H,psi;nsweeps=remaining,maxdim=cap,cutoff,noise=0,controls...,observer=obs,outputlevel=1)
    S,p=entropy_at(psi,lat.cut)
    result=Dict("family"=>family,"length"=>L,"width"=>w,"spins"=>lat.n,"vertices"=>lat.nv,"nup"=>nup,"explicit_number_sector"=>explicit_sector,
        "physical_circumference"=>lat.circumference,"ordering"=>ordering,"cut"=>lat.cut,"legs"=>lat.legs,
        "vertices_table"=>lat.verts,"mps_order"=>lat.order,"seed"=>seed,"audit"=>audit,"git_commit"=>audit["git_commit"],
        "runtime_seconds"=>time()-started,"resumed_from"=>path,"records"=>[Dict("cap"=>cap,"energy"=>E,"entropy"=>S,
        "maxlinkdim"=>maxlinkdim(psi),"schmidt_probabilities"=>p,"sweep_energies"=>energies(obs),"sweep_max_truncation_errors"=>truncerrors(obs))])
    manifest=completed_finite_checkpoint(psi,result,cap,"resumed")
    point="results/$(family)_L$(L)_w$(w)_chi$(cap)_$(ordering)_resumed.json"
    seed==7103 || (point=replace(point,".json"=>"_seed$(seed).json"))
    explicit_sector && (point=replace(point,".json"=>"_Nup$(nup).json"))
    atomic_json(point,manifest)
    current=manifest["checkpoint_file"]
    for nextcap in [32,64,128,256,512,1024,2048]
        cap<nextcap<=target || continue
        continue_cylinder(current,nextcap;cutoff=min(cutoff,1e-11),controls...)
        output="results/$(family)_L$(L)_w$(w)_chi$(nextcap)"*(ordering=="axial" ? "" : "_"*ordering)
        seed==7103 || (output*="_seed$(seed)")
        explicit_sector && (output*="_Nup$(nup)")
        current=JSON3.read(read(output*".json",String),Dict{String,Any})["checkpoint_file"]
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS[1]=="--latest"
        path=newest_checkpoint(ARGS[2],parse(Int,ARGS[3]),parse(Int,ARGS[4]))
        resume_checkpoint(path,length(ARGS)>=5 ? parse(Int,ARGS[5]) : nothing)
    else
        resume_checkpoint(ARGS[1],length(ARGS)>=2 ? parse(Int,ARGS[2]) : nothing)
    end
end
