# A logging/checkpoint driver around the package's own VUMPS iteration.
# Tensor contractions, environments, expansion, and eigensolves remain library code.
function audited_vumps(H,psi;family,w,cap,maxiter=30,tol=1e-7,
        solver_tol=x->1e-10,ordering="matter_first",tag="",seed=7103,audit=nothing)
    n=18w
    audit===nothing && (audit=run_provenance(;seed,solver="ITensorInfiniteMPS VUMPS",
        settings=Dict("tol"=>tol,"local_eigensolver_tolerance"=>"caller solver_tol","maxiter"=>maxiter),
        initialization="input canonical infinite MPS",conserved_quantum_numbers=hasqns(siteind(psi.AL,1)) ? ["total Sz"] : String[]))
    epsL=fill(tol,n);epsR=fill(tol,n)
    started=time();iterations=[]
    base="results/infinite_$(family)_w$(w)_chi$(cap)$(tag)"
    sector=hasqns(siteind(psi.AL,1)) ? "Sz0" : "unrestricted"
    checkpoint="results/checkpoints/infinite_$(family)_w$(w)_cell2_N$(n)_chi$(cap)_$(sector)_seed$(seed)_v$(RUN_FORMAT_VERSION)$(tag)_latest.jls"
    for iteration=1:maxiter
        elapsed=@elapsed psi,(left_energy,right_energy)=ITensorInfiniteMPS.tdvp_iteration(
            ITensorInfiniteMPS.vumps_solver,H,psi; (ϵᴸ!)=epsL,(ϵᴿ!)=epsR,
            multisite_update_alg="sequential",solver_tol,time_step=-Inf,eager=true)
        residual=max(maximum(epsL),maximum(epsR))
        record=Dict("iteration"=>iteration,"bond_dimension"=>maxlinkdim(psi),"canonical_solver_residual"=>residual,
            "runtime_seconds"=>elapsed,"left_energy_real"=>real.(left_energy),"left_energy_imag"=>imag.(left_energy),
            "right_energy_real"=>real.(right_energy),"right_energy_imag"=>imag.(right_energy),"tol"=>tol,
            "git_commit"=>audit["git_commit"],"run_version"=>RUN_FORMAT_VERSION)
        push!(iterations,record)
        open(base*"_iterations.jsonl","a") do io;println(io,JSON3.write(record));end
        println("AUDITED VUMPS iteration=",iteration," chi=",maxlinkdim(psi)," residual=",residual," seconds=",elapsed);flush(stdout)
        if iseven(iteration) || residual<tol || iteration==maxiter
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
        residual<tol && break
    end
    psi
end
