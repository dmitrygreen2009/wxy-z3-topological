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
    if get(metadata,"kind",nothing)=="finite" && haskey(metadata,"nup") && hasqns(psi)
        n=length(psi)
        @assert val(flux(psi),"Sz")==2metadata["nup"]-n "Checkpoint number sector does not match the physical MPS"
    end
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
            "basis"=>get(get(audit,"solver_settings",Dict()),"basis","physical_spin"),
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
    endswith(path,".previous") && isfile(path) && return path
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
function load_state(path)
    actual=resolve_checkpoint(path)
    sidecar=replace(actual,".jls"=>".json")
    if isfile(sidecar)
        metadata=JSON3.read(read(sidecar,String),Dict{String,Any})
        if haskey(metadata,"checkpoint_sha256") && haskey(metadata,"kind")
            return first(load_valid_checkpoint(actual))
        end
    end
    # Legacy payloads and the short finite-to-infinite fitting cache lack a
    # checksummed manifest; preserve their compatibility explicitly.
    deserialize(actual)
end
function completed_finite_checkpoint(psi,result,cap,phase)
    family=result["family"];L=result["length"];w=result["width"];n=result["spins"]
    seed=get(result,"seed",7103);nup=get(result,"nup",cld(n,2));Ly=round(result["physical_circumference"];digits=5)
    order=get(result,"ordering","axial")
    completion_id=replace(string(now(UTC)),r"[^0-9]"=>"")
    path="results/checkpoints/$(family)_w$(w)_Ly$(Ly)_L$(L)_N$(n)_chi$(cap)_Nup$(nup)_seed$(seed)_v$(RUN_FORMAT_VERSION)_$(order)_$(phase)_complete_$(completion_id).jls"
    basepath=path;collision=1
    while isfile(path)
        path=replace(basepath,".jls"=>"_$(collision).jls");collision+=1
    end
    metadata=merge(copy(result),Dict("kind"=>"finite","physical_spins"=>n,"cap"=>cap,"phase"=>phase))
    save_checkpoint(path,psi,metadata)
end
struct ResourceLimitError <: Exception
    message::String
end
Base.showerror(io::IO,e::ResourceLimitError)=print(io,e.message)
# The finite-library maxlinkdim omits the final unit-cell wrap bond.
# Inspect every periodic bond without changing tensors or solver behavior.
function mps_bond_dimension_details(psi)
    levels=hasproperty(psi,:AL) && hasproperty(psi,:AR) ? Dict("AL"=>psi.AL,"AR"=>psi.AR) : Dict("MPS"=>psi)
    periodic=hasproperty(psi,:AL) || (hasproperty(psi,:data) && hasproperty(getproperty(psi,:data),:translator))
    dimensions=Dict{String,Vector{Int}}()
    for (label,A) in levels
        count=periodic ? length(A) : max(0,length(A)-1)
        dimensions[label]=[let link=commonind(A[j],A[j+1]);isnothing(link) ? 1 : dim(link) end for j=1:count]
    end
    largest=maximum(vcat(collect(values(dimensions))...);init=1)
    Dict("periodic_wrap_bond_included"=>periodic,"bond_dimensions_by_tensor_representation"=>dimensions,
        "maximum_bond_dimension"=>largest,"historical_library_maxlinkdim"=>maxlinkdim(psi))
end
maximum_bond_dimension(psi)=mps_bond_dimension_details(psi)["maximum_bond_dimension"]

function estimate_memory(psi,target_chi;label="calculation")
    details=mps_bond_dimension_details(psi);current=max(1,details["maximum_bond_dimension"]);state_bytes=Base.summarysize(psi)
    # Conservative environment/Krylov allowance; record estimates separately from measurements.
    estimate=round(Int,1.5*2.0^30+20state_bytes*max(1.,target_chi/current)^2)
    budget=round(Int,0.5Sys.total_memory())
    record=Dict("label"=>label,"current_chi"=>current,"target_chi"=>target_chi,
        "bond_dimension_details"=>details,"state_bytes"=>state_bytes,"estimated_peak_bytes"=>estimate,"process_memory_budget_bytes"=>budget,
        "system_memory_bytes"=>Sys.total_memory(),"measured_process_peak_bytes"=>Sys.maxrss(),
        "policy"=>"1.5 GiB runtime plus 20 times projected MPS storage; one half of system RAM per large worker",
        "allowed"=>estimate<=budget,"git_commit"=>try readchomp(`git rev-parse HEAD`) catch; "unknown" end)
    open("results/resource_estimates.jsonl","a") do io;println(io,JSON3.write(record));end
    estimate<=budget || throw(ResourceLimitError("Predicted peak $(round(estimate/2.0^30;digits=2)) GiB exceeds the $(round(budget/2.0^30;digits=2)) GiB per-worker budget for $label; existing checkpoints preserved."))
    record
end
