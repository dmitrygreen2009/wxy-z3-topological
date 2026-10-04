include("cylinders.jl")
include("continue.jl")
include("infinite.jl")
function newest_checkpoint(family,L,w)
    candidates=[]
    for file in readdir("results/checkpoints";join=true)
        endswith(file,".json") || endswith(file,".json.previous") || continue
        meta=try JSON3.read(read(file,String),Dict{String,Any}) catch; continue end
        get(meta,"family",nothing)==family && get(meta,"length",nothing)==L && get(meta,"width",nothing)==w || continue
        get(meta,"kind",nothing)=="finite" || continue
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
        H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering),siteinds(only,psi))
        psi=audited_vumps(H,psi;family,w,cap,maxiter=max(1,get(meta,"stage_maxiter",40)-meta["iteration"]),ordering,tag)
        measure_infinite(psi,H,family,w,cap,meta["iteration"]+1;ordering,tag)
        return
    end
    L=meta["length"];seed=get(meta,"random_seed",7103);ordering=meta["ordering"]
    lat=cylinder(family,L,w;ordering);H=MPO(lat.os,siteinds(psi))
    phase=string(get(meta,"stage",get(meta,"phase","resume")))
    planned=occursin("refinement",phase) ? 6 : 8
    remaining=max(1,planned-get(meta,"sweep",0))
    target===nothing && (target=get(get(meta,"solver_settings",Dict()),"target_chi",cap))
    audit=run_provenance(;seed,solver="ITensorMPS DMRG resume",settings=Dict("remaining_sweeps"=>remaining,"cutoff"=>get(meta,"cutoff",1e-11),"target_chi"=>target),
        initialization=path,conserved_quantum_numbers=["total Sz"])
    cutoff=get(meta,"cutoff",1e-11)
    obs=finite_checkpoint_observer(lat;family,L,w,cap,seed,stage="resumed",ordering,audit,cutoff,noise=0)
    started=time()
    E,psi=dmrg(H,psi;nsweeps=remaining,maxdim=cap,cutoff,noise=0,eigsolve_krylovdim=12,observer=obs,outputlevel=1)
    S,p=entropy_at(psi,lat.cut)
    result=Dict("family"=>family,"length"=>L,"width"=>w,"spins"=>lat.n,"vertices"=>lat.nv,
        "physical_circumference"=>lat.circumference,"ordering"=>ordering,"cut"=>lat.cut,"legs"=>lat.legs,
        "vertices_table"=>lat.verts,"mps_order"=>lat.order,"seed"=>seed,"audit"=>audit,"git_commit"=>audit["git_commit"],
        "runtime_seconds"=>time()-started,"resumed_from"=>path,"records"=>[Dict("cap"=>cap,"energy"=>E,"entropy"=>S,
        "maxlinkdim"=>maxlinkdim(psi),"schmidt_probabilities"=>p,"sweep_energies"=>energies(obs),"sweep_max_truncation_errors"=>truncerrors(obs))])
    manifest=completed_finite_checkpoint(psi,result,cap,"resumed")
    point="results/$(family)_L$(L)_w$(w)_chi$(cap)_$(ordering)_resumed.json"
    atomic_json(point,manifest)
    current=manifest["checkpoint_file"]
    for nextcap in [32,64,128,256,512,1024,2048]
        cap<nextcap<=target || continue
        continue_cylinder(current,nextcap)
        output="results/$(family)_L$(L)_w$(w)_chi$(nextcap)"*(ordering=="axial" ? "" : "_"*ordering)
        seed==7103 || (output*="_seed$(seed)")
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
