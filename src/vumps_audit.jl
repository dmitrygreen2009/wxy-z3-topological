include("infinite_number_background.jl")
include("local_eigenpair_audit.jl")
# A logging/checkpoint driver around the package's own VUMPS iteration.
# Tensor contractions, environments, expansion, and eigensolves remain library code.
function audited_vumps(H,psi;family,w,cap,maxiter=30,tol=1e-7,
        solver_tol=x->1e-10,solver_tolerance_rule="fixed",ordering="matter_first",tag="",seed=7103,audit=nothing,update_algorithm="sequential")
    @assert update_algorithm in ["sequential","parallel"]
    n=18w
    audit===nothing && (audit=run_provenance(;seed,solver="ITensorInfiniteMPS VUMPS",
        settings=Dict("tol"=>tol,"local_eigensolver_tolerance_at_initial_residual"=>solver_tol(tol),"local_solver_tolerance_rule"=>solver_tolerance_rule,"multisite_update_algorithm"=>"sequential","time_step"=>"-Inf","subspace_expansion_cutoff"=>1e-10,"maxiter"=>maxiter),
        initialization="input canonical infinite MPS",conserved_quantum_numbers=hasqns(siteind(psi.AL,1)) ? ["total Sz"] : String[]))
    epsL=fill(tol,n);epsR=fill(tol,n)
    started=time();iterations=[]
    settings=Dict{String,Any}(get(audit,"solver_settings",Dict()))
    settings["multisite_update_algorithm"]=update_algorithm;audit["solver_settings"]=settings
    base="results/infinite_$(family)_w$(w)_chi$(cap)$(tag)"
    number_background=infinite_number_background(psi);audit["number_background"]=number_background
    sector=number_background["label"]
    checkpoint="results/checkpoints/infinite_$(family)_w$(w)_cell2_N$(n)_chi$(cap)_$(sector)_seed$(seed)_v$(RUN_FORMAT_VERSION)$(tag)_latest.jls"
    audit["local_solver_evidence"]="Library ConvergenceInfo and independently evaluated normalized eigenpair residuals; prior records without these fields remain unaudited locally."
    for iteration=1:maxiter
        requested_local_tolerance=solver_tol(max(maximum(epsL),maximum(epsR)))
        local_records=[]
        function logged_solver(M,time_step,v0,local_tol,eager=true)
            lambda,v,info=ITensorInfiniteMPS.vumps_solver(M,time_step,v0,local_tol,eager)
            push!(local_records,local_eigenpair_record(M,lambda,v,info,local_tol))
            lambda,v,info
        end
        elapsed=@elapsed psi,(left_energy,right_energy)=ITensorInfiniteMPS.tdvp_iteration(
            logged_solver,H,psi; (ϵᴸ!)=epsL,(ϵᴿ!)=epsR,
            multisite_update_alg=update_algorithm,solver_tol,time_step=-Inf,eager=true)
        residual=max(maximum(epsL),maximum(epsR))
        local_passed=all(r["local_criterion_passed"] for r in local_records)
        open(base*"_local_eigensolves.jsonl","a") do io
            println(io,JSON3.write(Dict("iteration"=>iteration,"git_commit"=>audit["git_commit"],"solves"=>local_records)))
        end
        record=Dict("iteration"=>iteration,"bond_dimension"=>maxlinkdim(psi),"canonical_solver_residual"=>residual,
            "runtime_seconds"=>elapsed,"left_energy_real"=>real.(left_energy),"left_energy_imag"=>imag.(left_energy),
            "right_energy_real"=>real.(right_energy),"right_energy_imag"=>imag.(right_energy),"tol"=>tol,
            "local_eigensolver_tolerance_evaluated"=>requested_local_tolerance,"git_commit"=>audit["git_commit"],"run_version"=>RUN_FORMAT_VERSION,
            "local_solve_count"=>length(local_records),"all_local_eigenpairs_converged"=>local_passed,
            "max_true_local_eigenpair_residual"=>maximum(r["true_eigenpair_residual"] for r in local_records),
            "minimum_local_converged_eigenpairs"=>minimum(r["converged_eigenpairs"] for r in local_records))
        push!(iterations,record)
        open(base*"_iterations.jsonl","a") do io;println(io,JSON3.write(record));end
        println("AUDITED VUMPS iteration=",iteration," chi=",maxlinkdim(psi)," residual=",residual," seconds=",elapsed);flush(stdout)
        if iseven(iteration) || (residual<tol && local_passed) || iteration==maxiter
            metadata=merge(copy(audit),Dict("kind"=>"infinite","family"=>family,"width"=>w,
                "physical_circumference"=>family=="zigzag" ? sqrt(3)*w : 3.0w,
                "cell_slices"=>2,"cell_spins"=>n,"cap"=>cap,"bond_dimension"=>maxlinkdim(psi),
                "infinite_ordering"=>ordering,"measurement_tag"=>tag,"iteration"=>iteration,
                "stage_maxiter"=>maxiter,"tol"=>tol,"canonical_solver_residual"=>residual,
                "runtime_seconds"=>time()-started,"iterations"=>iterations,
                "checkpoint_interpretation"=>"In-progress variational state; perform full transfer and entropy checks before scientific use."))
            save_checkpoint(checkpoint,psi,metadata)
            atomic_json(base*"_solver_audit.json",metadata)
        end
        residual<tol && local_passed && break
    end
    psi
end
