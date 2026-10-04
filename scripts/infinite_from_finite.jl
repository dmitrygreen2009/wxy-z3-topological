include("infinite.jl")
function infinite_from_finite(path;initial_chi=16,target_chi=64)
    Random.seed!(7103)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    family=meta["family"];w=meta["width"];L=meta["length"]
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
    if isfile(fit_path)
        psi=load_state(fit_path)
    else
        Random.seed!(7103)
        psi=infinitemps_approx(finite;nsites=n,nrange=middle,nsweeps=4,outputlevel=1)
        serialize(fit_path,psi)
    end
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
    for (iteration,cap) in enumerate(unique([initial_chi,min(32,target_chi),target_chi]))
        estimate_memory(psi,cap;label="$(family) infinite w$(w) expansion")
        psi=subspace_expansion(psi,H;cutoff=1e-10,maxdim=cap)
        psi=audited_vumps(H,psi;family,w,cap,tag="_finite_seed_star",ordering="star",tol=1e-7,maxiter=40,
            solver_tol=x->1e-10)
        measure_infinite(psi,H,family,w,cap,iteration;tag="_finite_seed_star",ordering="star")
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;infinite_from_finite(ARGS[1]);end
