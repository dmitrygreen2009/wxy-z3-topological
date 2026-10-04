include("infinite.jl")
include("../src/infinite_block_initializer.jl")
function run_random_cell(family,w,nup;initial_chi=8,target_chi=64,seed=7226,cell_slices=2,update_algorithm="sequential",ordering="star")
    n=9w*cell_slices;@assert 0<nup<n
    @assert family in ["zigzag","armchair"]
    initstate(j)=fld(mod1(j,n)*nup,n)>fld((mod1(j,n)-1)*nup,n) ? "Up" : "Dn"
    @assert count(j->initstate(j)=="Up",1:n)==nup
    Random.seed!(seed)
    s=infsiteinds("S=1/2",n;conserve_qns=true,initstate)
    phi=random_mps(ComplexF64,collect(s[Cell(1)]),[initstate(j) for j=1:n];linkdims=initial_chi)
    @assert iszero(flux(phi))
    psi=infinite_cell_product(phi;tol=1e-14)
    density=nup//n
    tag="_qn_rho$(numerator(density))of$(denominator(density))_cell_product_seed$(seed)_$(ordering)"
    cell_slices!=2 && (tag*="_cell$(cell_slices)")
    initial_path="results/checkpoints/infinite_$(family)_w$(w)_cell$(cell_slices)_N$(n)_chi$(initial_chi)_rho$(numerator(density))of$(denominator(density))_seed$(seed)_v$(RUN_FORMAT_VERSION)$(tag)_initial.jls"
    if isfile(initial_path)
        psi=load_state(initial_path)
    else
        save_checkpoint(initial_path,psi,merge(run_provenance(seed=seed,solver="ITensorInfiniteMPS independent-cell initializer",
            initialization="Fixed-density library random_mps with deliberate unit-dimensional cell boundary",
            conserved_quantum_numbers=["U1 mean physical density $(density)"]),
            Dict("kind"=>"infinite_initializer","family"=>family,"width"=>w,"cell_spins"=>n,"cell_slices"=>cell_slices,"infinite_ordering"=>ordering,"measurement_tag"=>tag,
                "number_background"=>infinite_number_background(psi),"checkpoint_interpretation"=>"Initializer only; not an optimized scientific result")))
    end
    H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering,cell_slices),siteinds(only,psi))
    caps=sort(unique(vcat(initial_chi,[c for c in [16,32,64,128,256] if initial_chi<c<target_chi],target_chi)))
    for (stage,cap) in enumerate(caps)
        stage_seed=seed+stage;Random.seed!(stage_seed)
        estimate_memory(psi,cap;label="$(family) w$(w) fixed-density random-cell QN initialization")
        psi=subspace_expansion(psi,H;cutoff=1e-10,maxdim=cap)
        audit=run_provenance(seed=stage_seed,solver="ITensorInfiniteMPS QN VUMPS from independent entangled-cell initializer",
            settings=Dict("tol"=>1e-7,"maxiter"=>40,"local_eigensolver_tolerance"=>1e-10,
                "initial_chi"=>initial_chi,"subspace_expansion_cutoff"=>1e-10,"canonicalization_tolerance"=>1e-14,
                "target_physical_nup_per_cell"=>nup,"target_physical_number_density"=>nup/n,
                "individual_physical_sites_per_cell"=>n,"initial_intercell_virtual_dimension"=>1,
                "physical_site_QN_spaces"=>[string(space(j)) for j in s[Cell(1)]],"initial_state_seed"=>seed),
            initialization="Normalized library random_mps in a fixed-number finite cell, repeated with deliberate unit-dimensional cell boundaries; no Hamiltonian approximation",
            conserved_quantum_numbers=["U(1), physical number density $nup/$n via official shifted/scaled site QNs"])
        psi=audited_vumps(H,psi;family,w,cap,tag,ordering,tol=1e-7,maxiter=40,solver_tol=x->1e-10,audit,seed=stage_seed,update_algorithm)
        measure_infinite(psi,H,family,w,cap,stage;tag,ordering)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    run_random_cell(ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3]);target_chi=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 64,cell_slices=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 2,update_algorithm=length(ARGS)>=6 ? ARGS[6] : "sequential",ordering=length(ARGS)>=7 ? ARGS[7] : "star")
end
