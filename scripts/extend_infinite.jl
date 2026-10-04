include("infinite.jl")
for family in ["zigzag","armchair"]
    stem="results/infinite_$(family)_w1_chi64"
    psi=load_state(stem*".jls")
    ss=siteinds(only,psi);H=InfiniteSum{MPO}(infinite_opsum(family,1),ss)
    psi=subspace_expansion(psi,H;cutoff=1e-9,maxdim=128)
    psi=tdvp(H,psi;time_step=-Inf,tol=1e-7,maxiter=60,
        solver_tol=x->max(x/10,1e-10),multisite_update_alg="sequential")
    measure_infinite(psi,H,family,1,128,8)
    run_infinite(family,2,64)
end
