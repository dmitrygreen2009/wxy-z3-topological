using Dates, SHA, Serialization
const RUN_FORMAT_VERSION=3
# Capture once when this process loads the project, not after files are edited
# while a long calculation is running.
const LAUNCH_REVISION=try readchomp(`git rev-parse HEAD`) catch; "unknown" end
const LAUNCH_CODE_DIRTY=try !isempty(readchomp(`git status --porcelain -- src scripts test Project.toml Manifest.toml requirements.txt`)) catch; true end
const LAUNCH_SOURCE_HASHES=let fingerprints=Dict()
    for root in ["src","scripts","test"],file in readdir(root;join=true)
        isfile(file) || continue
        fingerprints[file]=bytes2hex(sha256(read(file)))
    end
    fingerprints
end
const PROCESS_STARTED_UTC=string(now(UTC))
function run_provenance(;seed=nothing,solver,settings=Dict(),initialization,conserved_quantum_numbers)
    revision=LAUNCH_REVISION;dirty=LAUNCH_CODE_DIRTY;fingerprints=copy(LAUNCH_SOURCE_HASHES)
    Dict("run_version"=>RUN_FORMAT_VERSION,"git_commit"=>revision,"code_worktree_dirty"=>dirty,
        "source_sha256"=>fingerprints,"process_started_utc"=>PROCESS_STARTED_UTC,"provenance_snapshot"=>"process launch","julia_version"=>string(VERSION),"started_utc"=>string(now(UTC)),
        "random_seed"=>seed,"solver"=>solver,"solver_settings"=>settings,"initialization"=>initialization,
        "conserved_quantum_numbers"=>conserved_quantum_numbers)
end
function atomic_json(path,data)
    mkpath(dirname(path));temporary=path*".tmp"
    open(temporary,"w") do io;JSON3.write(io,data);end
    mv(temporary,path;force=true)
end
function save_checkpoint(path,psi,metadata)
    mkpath(dirname(path));temporary=path*".tmp"
    serialize(temporary,psi)
    if isfile(path)
        mv(path,path*".previous";force=true)
        sidecar=replace(path,".jls"=>".json")
        isfile(sidecar) && mv(sidecar,sidecar*".previous";force=true)
    end
    mv(temporary,path;force=true)
    manifest=copy(metadata)
    manifest["checkpoint_file"]=path;manifest["checkpoint_bytes"]=filesize(path)
    manifest["checkpoint_written_utc"]=string(now(UTC))
    manifest["checkpoint_sha256"]=open(io->bytes2hex(sha256(io)),path)
    atomic_json(replace(path,".jls"=>".json"),manifest)
    manifest
end
function load_valid_checkpoint(path)
    metadata=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    @assert filesize(path)==metadata["checkpoint_bytes"]
    @assert open(io->bytes2hex(sha256(io)),path)==metadata["checkpoint_sha256"]
    psi=deserialize(path)
    if metadata["kind"]=="finite"
        @assert isfinite(norm(psi)) && abs(norm(psi)-1)<1e-6
    else
        @assert all(isfinite(norm(psi.C[j])) && abs(norm(psi.C[j])-1)<1e-6 for j=1:metadata["cell_spins"])
    end
    psi,metadata
end
function finite_checkpoint_observer(lat;family,L,w,cap,seed,stage,ordering,audit,cutoff,noise,nup=cld(lat.n,2))
    started=time();holder=Ref{Any}()
    stem="results/checkpoints/$(family)_w$(w)_L$(L)_N$(lat.n)_chi$(cap)_Nup$(nup)_seed$(seed)_v$(RUN_FORMAT_VERSION)_$(ordering)_$(stage)_latest.jls"
    callback=function(kwargs)
        sweep=kwargs[:sweep]
        iseven(sweep) || return
        psi=kwargs[:psi]
        # Entropy measurement must not alter the live DMRG state or environments.
        S,p=entropy_at(deepcopy(psi),lat.cut)
        metadata=merge(copy(audit),Dict("kind"=>"finite","family"=>family,"length"=>L,"width"=>w,
            "physical_circumference"=>lat.circumference,"physical_spins"=>lat.n,"ordering"=>ordering,
            "cap"=>cap,"bond_dimension"=>maxlinkdim(psi),"nup"=>nup,"stage"=>stage,"sweep"=>sweep,
            "energy"=>real(kwargs[:energy]),"full_state_entropy"=>S,"schmidt_probabilities"=>p,
            "cutoff"=>cutoff,"noise"=>noise,"runtime_seconds"=>time()-started,
            "sweep_energies"=>energies(holder[]),"sweep_max_truncation_errors"=>truncerrors(holder[]),
            "correlation_length"=>nothing,"correlation_length_status"=>"Finite MPS; no arbitrary periodic bond identification"))
        save_checkpoint(stem,psi,metadata)
    end
    observer=LoggedObserver(;checkpoint_callback=callback);holder[]=observer
    observer
end

function resolve_checkpoint(path)
    sidecar=replace(path,".jls"=>".json")
    if isfile(sidecar)
        metadata=JSON3.read(read(sidecar,String),Dict{String,Any})
        if haskey(metadata,"checkpoint_file") && isfile(metadata["checkpoint_file"])
            return metadata["checkpoint_file"]
        end
    end
    @assert isfile(path) "No valid checkpoint or result manifest for $path"
    path
end
load_state(path)=deserialize(resolve_checkpoint(path))
function completed_finite_checkpoint(psi,result,cap,phase)
    family=result["family"];L=result["length"];w=result["width"];n=result["spins"]
    seed=get(result,"seed",7103);nup=get(result,"nup",cld(n,2));Ly=round(result["physical_circumference"];digits=5)
    order=get(result,"ordering","axial")
    path="results/checkpoints/$(family)_w$(w)_Ly$(Ly)_L$(L)_N$(n)_chi$(cap)_Nup$(nup)_seed$(seed)_v$(RUN_FORMAT_VERSION)_$(order)_$(phase)_complete.jls"
    metadata=merge(copy(result),Dict("kind"=>"finite","physical_spins"=>n,"cap"=>cap,"phase"=>phase))
    save_checkpoint(path,psi,metadata)
end
struct ResourceLimitError <: Exception
    message::String
end
Base.showerror(io::IO,e::ResourceLimitError)=print(io,e.message)
function estimate_memory(psi,target_chi;label="calculation")
    current=max(1,maxlinkdim(psi));state_bytes=Base.summarysize(psi)
    # Conservative environment/Krylov allowance; record estimates separately from measurements.
    estimate=round(Int,1.5*2.0^30+20state_bytes*(target_chi/current)^2)
    budget=round(Int,0.5Sys.total_memory())
    record=Dict("label"=>label,"current_chi"=>current,"target_chi"=>target_chi,
        "state_bytes"=>state_bytes,"estimated_peak_bytes"=>estimate,"process_memory_budget_bytes"=>budget,
        "system_memory_bytes"=>Sys.total_memory(),"measured_process_peak_bytes"=>Sys.maxrss(),
        "policy"=>"1.5 GiB runtime plus 20 times projected MPS storage; one half of system RAM per large worker",
        "allowed"=>estimate<=budget,"git_commit"=>try readchomp(`git rev-parse HEAD`) catch; "unknown" end)
    open("results/resource_estimates.jsonl","a") do io;println(io,JSON3.write(record));end
    estimate<=budget || throw(ResourceLimitError("Predicted peak $(round(estimate/2.0^30;digits=2)) GiB exceeds the $(round(budget/2.0^30;digits=2)) GiB per-worker budget for $label; existing checkpoints preserved."))
    record
end
