# Read/setup only: never calls DMRG or a VUMPS iteration.
include("resume_checkpoint.jl")
plans=JSON3.read(read(ARGS[1],String),Dict{String,Any})
records=[]
for job in plans["jobs"]
    started=time();record=Dict{String,Any}("pid"=>job["pid"],"checkpoint"=>job["snapshot"],"purpose"=>job["purpose"])
    try
        record["validation"]=resume_checkpoint(job["snapshot"],job["target_cap"];validate_only=true,preserve_stage_budget=true)
        record["validated"]=true
    catch error
        record["validated"]=false;record["error"]=sprint(showerror,error,catch_backtrace())
    end
    record["runtime_seconds"]=time()-started
    push!(records,record)
    atomic_json("results/julia_checkpoint_restart_validation.json",Dict("records"=>records,"julia_version"=>string(VERSION),"no_optimization_executed"=>true))
    println("Checkpoint resume setup ",job["pid"],": ",record["validated"]);flush(stdout)
    GC.gc()
end
# Read the finite seed separately, without claiming it reconstructs its live batch.
seed=plans["protected_finite_job"]["checkpoint_payload"]
try
    psi=load_state(seed)
    atomic_json("results/julia_finite_seed_readability.json",Dict("checkpoint"=>seed,"readable"=>true,"spins"=>length(psi),
        "maxlinkdim"=>maxlinkdim(psi),"protected_process"=>plans["protected_finite_job"]["pid"],
        "complete_batch_resume_validated"=>false,"decision"=>"Keep live process paused: current higher-bond sweep and batch control state are not in this seed"))
catch error
    atomic_json("results/julia_finite_seed_readability.json",Dict("checkpoint"=>seed,"readable"=>false,"error"=>sprint(showerror,error),"complete_batch_resume_validated"=>false))
end
