include("infinite.jl")
function infinite_from_finite(path;initial_chi=16,target_chi=64,tag="_finite_seed_star")
    Random.seed!(7103)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    family=meta["family"];w=meta["width"];L=meta["length"]
    @assert L>=4 && iseven(L) "Bulk two-slice warm-start indexing requires an even length of at least four"
    @assert get(meta,"ordering","axial")=="star"
    lat=cylinder(family,L,w;ordering="star");n=18w
    finite=dense(load_state(path))
    truncate!(finite;maxdim=initial_chi,cutoff=1e-10)
    # Match exactly two bulk spatial slices to the infinite star ordering.
    leading=lat.n-9w*L
    t=clamp(fld(L-2,2),1,L-3)
    middle=leading+9w*t+1:leading+9w*t+n
    @assert length(middle)==n
    println("INFINITE FIT physical bulk slices ",t,":",t+1," MPS range ",middle);flush(stdout)
    # Official library variational fitting; unrelated finite bonds are not identified.
    fit_path=replace(path,".jls"=>"_infinite_fit_chi$(initial_chi).jls")
    source_payload=resolve_checkpoint(path)
    source_sha=open(io->bytes2hex(sha256(io)),source_payload)
    cache_metadata_path=replace(fit_path,".jls"=>"_cache_audit.json")
    cache_metadata=isfile(cache_metadata_path) ? JSON3.read(read(cache_metadata_path,String),Dict{String,Any}) : Dict{String,Any}()
    if isfile(fit_path) && get(cache_metadata,"fit_source_checkpoint_sha256",nothing)!==nothing && cache_metadata["fit_source_checkpoint_sha256"]!=source_sha
        # This cache is demonstrably from another physical input. Preserve it,
        # and use a distinct cache for the changed source instead of overwriting.
        fit_path=replace(fit_path,".jls"=>"_source$(source_sha[1:12]).jls")
        cache_metadata_path=replace(fit_path,".jls"=>"_cache_audit.json")
        cache_metadata=isfile(cache_metadata_path) ? JSON3.read(read(cache_metadata_path,String),Dict{String,Any}) : Dict{String,Any}()
    end
    @assert get(cache_metadata,"fit_source_checkpoint_sha256",source_sha) in (nothing,source_sha) "Fitting cache input fingerprint mismatch"
    cached=isfile(fit_path)
    if cached
        if haskey(cache_metadata,"fit_payload_sha256")
            @assert open(io->bytes2hex(sha256(io)),fit_path)==cache_metadata["fit_payload_sha256"] "Fitting cache checksum mismatch; preserve the payload and diagnose it"
        end
        psi=load_state(fit_path)
    else
        Random.seed!(7103)
        psi=infinitemps_approx(finite;nsites=n,nrange=middle,nsweeps=4,outputlevel=1)
        serialize(fit_path*".tmp",psi);mv(fit_path*".tmp",fit_path;force=true)
        cache_metadata["fit_source_checkpoint_sha256"]=source_sha
    end
    cache_metadata["requested_finite_result"]=path
    cache_metadata["requested_source_checkpoint_sha256"]=source_sha
    cache_metadata["fit_payload_sha256"]=open(io->bytes2hex(sha256(io)),fit_path)
    cache_metadata["fit_payload_bytes"]=filesize(fit_path)
    cache_metadata["fit_source_status"]=get(cache_metadata,"fit_source_checkpoint_sha256",nothing)===nothing ? "Legacy origin unverified; reusable trial initialization, not claimed to originate from the latest requested finite state" : "Source fingerprint recorded"
    cache_metadata["requested_fit_settings"]=Dict("initial_chi"=>initial_chi,"finite_truncation_cutoff"=>1e-10,"fit_sweeps"=>4,"random_seed"=>7103,"bulk_range"=>collect(middle))
    cache_metadata["recorded_git_commit"]=LAUNCH_REVISION
    atomic_json(cache_metadata_path,cache_metadata)
    # mixed_canonical's current wrapper does not forward its tol keyword.
    # Call the official right/left canonicalization routines with stricter tol.
    trials=Dict{String,Any}[]
    right=nothing
    # The package rejects a positive imaginary eigenvalue component >1e-15,
    # even when the eigensolver converges at 1e-14. Retry its random start;
    # preserve the fitted state and every tolerance, and record each attempt.
    for attempt=1:8
        canonical_seed=7102+attempt;Random.seed!(canonical_seed)
        try
            _,right,_=ITensorInfiniteMPS.right_orthogonalize(psi.AL;left_tags=ts"Left",right_tags=ts"Right",tol=1e-14)
            push!(trials,Dict("seed"=>canonical_seed,"status"=>"success"));break
        catch error
            message=sprint(showerror,error)
            push!(trials,Dict("seed"=>canonical_seed,"status"=>"failure","error"=>message))
            atomic_json(replace(fit_path,".jls"=>"_canonicalization.json"),Dict("trials"=>trials,"tolerance"=>1e-14))
            occursin("Imaginary part of eigenvalue is large",message) || rethrow()
            attempt==8 && rethrow()
            println("Retry official canonicalization with a new seeded Krylov start: ",message);flush(stdout)
        end
    end
    atomic_json(replace(fit_path,".jls"=>"_canonicalization.json"),Dict("trials"=>trials,"tolerance"=>1e-14))
    left,center,lambda=ITensorInfiniteMPS.left_orthogonalize(right;tol=1e-14)
    @assert abs(lambda-1)<1e-10
    psi=InfiniteCanonicalMPS(left,center,right)
    ss=siteinds(only,psi);H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering="star"),ss)
    @assert 1<=initial_chi<=target_chi
    caps=sort(unique(vcat(initial_chi,[c for c in [16,32,64,128,256,512] if initial_chi<c<target_chi],target_chi)))
    for (iteration,cap) in enumerate(caps)
        stage_seed=7103+iteration;Random.seed!(stage_seed)
        estimate_memory(psi,cap;label="$(family) infinite w$(w) expansion")
        psi=subspace_expansion(psi,H;cutoff=1e-10,maxdim=cap)
        stage_audit=run_provenance(;seed=stage_seed,solver="ITensorInfiniteMPS VUMPS from finite-state fit",
            settings=Dict("tol"=>1e-7,"maxiter"=>40,"local_solver_tolerance_rule"=>"fixed","local_eigensolver_tolerance_at_initial_residual"=>1e-10,
                "time_step"=>"-Inf","multisite_update_algorithm"=>"sequential","subspace_expansion_cutoff"=>1e-10),
            initialization="Official infinite-state fit cache $fit_path; requested finite input $path",conserved_quantum_numbers=String[])
        stage_audit["warm_start_cache_audit"]=cache_metadata
        stage_audit["warm_start_source_projection_records"]=get(meta,"projection_records",nothing)
        stage_audit["infinite_CGS_charge_enforced"]=false
        psi=audited_vumps(H,psi;family,w,cap,tag,ordering="star",tol=1e-7,maxiter=40,audit=stage_audit,seed=stage_seed,
            solver_tol=x->1e-10)
        measure_infinite(psi,H,family,w,cap,iteration;tag,ordering="star")
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;infinite_from_finite(ARGS[1]);end
