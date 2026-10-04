# Exact commuting-projector penalty, never a gauge-only approximation.
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/cgs.jl")
include("../src/cgs_sector_penalty.jl")
function run_sector_penalty(family,L,w,nup,charge;seed=7252,maxcap=256,resume=nothing)
    started=time();ordering="star";Random.seed!(seed)
    lat=cylinder(family,L,w;ordering)
    @assert 0<=nup<=lat.n
    table="geometry/cgs_cycles/$(family)_L$(L)_w$(w)_star.json"
    definitions=JSON3.read(read(table,String),Dict{String,Any})
    geometry=JSON3.read(read(definitions["parent_geometry"],String),Dict{String,Any})
    @assert geometry["physical_spins"]==lat.n
    @assert Int.(geometry["mps_order_physical_spin_ids"])==lat.order
    @assert geometry["spatial_cut_mps_bond"]==lat.cut
    cycle=first(c for c in definitions["cycles"] if c["winding_number"]==1)
    p=Int.(cycle["edge_exponents"])
    if resume===nothing
        sites=siteinds("S=1/2",lat.n;conserve_qns=true)
        state=[j<=nup ? "Up" : "Dn" for j=1:lat.n]
        psi=random_mps(ComplexF64,sites,state;linkdims=4)
        psi,_=cgs_project(psi,lat,p;charge,maxdim=maxcap,cutoff=1e-14)
    else
        psi=load_state(resume);sites=siteinds(psi)
        @assert length(psi)==lat.n && val(flux(psi),"Sz")==2nup-lat.n
        previous=JSON3.read(read(replace(resume,".jls"=>".json"),String),Dict{String,Any})
        settings=get(previous,"solver_settings",get(get(previous,"audit",Dict()),"solver_settings",Dict()))
        @assert settings["winding_charge"]==charge && Int.(settings["edge_exponents"])==p
    end
    H=MPO(lat.os,sites);penalty=cgs_sector_penalty(lat,sites,p;charge)
    audit=run_provenance(;seed,solver="ITensorMPS DMRG sum of exact MPOs: H + lambda(I-P_q)",
        settings=Dict("winding_charge"=>charge,"edge_exponents"=>p,"penalty_strength"=>penalty.strength,
            "Hamiltonian_norm_bound"=>penalty.hamiltonian_norm_bound,"cutoff"=>1e-13,
            "sweeps_per_stage"=>12,"krylovdim"=>32,"eigsolve_maxiter"=>10,"eigsolve_tol"=>1e-13,
            "purity_error_tolerance"=>1e-10,"energy_variance_tolerance"=>1e-8,
            "energy_drift_tolerance"=>1e-9,"entropy_drift_tolerance"=>1e-5),
        initialization=resume===nothing ? "Random fixed-number MPS with CGS projector initializer" : resume,
        conserved_quantum_numbers=["N_up=$nup","target winding eigenvalue omega^$charge via exact commuting penalty"])
    stem="results/$(family)_L$(L)_w$(w)_Nup$(nup)_winding$(charge)_penalty_seed$(seed)"
    records=[];previousS=nothing;passed=false
    caps=unique(min.(maxcap,[32,64,128,256,512]))
    for cap in caps,pass in 1:2
        estimate_memory(psi,cap;label="Exact winding penalty $family L$L w$w q$charge")
        obs=finite_checkpoint_observer(lat;family,L,w,cap,seed,stage="winding$(charge)_penalty_pass$(pass)",ordering,audit,cutoff=1e-13,noise=0,nup)
        Epen,psi=dmrg(vcat([H],penalty.terms),psi;nsweeps=12,maxdim=cap,cutoff=1e-13,noise=0,
            eigsolve_krylovdim=32,eigsolve_maxiter=10,eigsolve_tol=1e-13,observer=obs,outputlevel=1)
        E=real(inner(psi',H,psi)/inner(psi,psi))
        mean=inner(psi',penalty.U,psi)/inner(psi,psi)
        variance=real(inner(H,psi,H,psi)/inner(psi,psi))-E^2
        @assert variance>=-1e-8 "Energy variance is numerically invalid"
        purity=abs(mean-cis(2pi*charge/3));S,probs=entropy_at(psi,lat.cut)
        energies_stage=energies(obs);drift=length(energies_stage)>1 ? abs(energies_stage[end]-energies_stage[end-1]) : Inf
        discarded=maximum(truncerrors(obs))
        entropy_drift=previousS===nothing ? nothing : abs(S-previousS)
        passed=purity<1e-10 && variance<1e-8 && drift<1e-9 && discarded<1e-10 && entropy_drift!==nothing && entropy_drift<1e-5
        @assert val(flux(psi),"Sz")==2nup-lat.n
        push!(records,Dict("cap"=>cap,"pass"=>pass,"energy"=>E,"penalized_energy"=>Epen,"entropy"=>S,
            "maxlinkdim"=>maxlinkdim(psi),"schmidt_probabilities"=>probs,"sweep_energies"=>energies_stage,
            "sweep_max_truncation_errors"=>truncerrors(obs),"energy_variance"=>variance,
            "loop_expectation_real"=>real(mean),"loop_expectation_imag"=>imag(mean),"loop_variance"=>1-abs2(mean),
            "sector_purity_error"=>purity,"entropy_change_previous_stage"=>entropy_drift,"converged_within_fixed_number_and_loop_sector"=>passed))
        result=Dict("family"=>family,"length"=>L,"width"=>w,"physical_circumference"=>lat.circumference,
            "spins"=>lat.n,"vertices"=>lat.nv,"ordering"=>ordering,"seed"=>seed,"nup"=>nup,"total_sz"=>nup-lat.n/2,
            "filling_fraction"=>nup/lat.n,"filling_selection"=>"Fixed conserved QN; not global filling minimization",
            "cut"=>lat.cut,"geometry_cycle_definitions"=>table,"loop_label"=>cycle["label"],"records"=>records,
            "physical_geometry"=>definitions["parent_geometry"],"cycle_definitions_sha256"=>bytes2hex(sha256(read(table))),
            "audit"=>audit,"git_commit"=>audit["git_commit"],"runtime_seconds"=>time()-started,
            "correlation_length"=>nothing,"correlation_length_status"=>"Finite MPS; no arbitrary infinite bond identification",
            "interpretation"=>"Ground candidate within the recorded number and exact microscopic winding sector. No MES, global filling or two-dimensional phase certificate.")
        atomic_json(stem*".json",completed_finite_checkpoint(psi,result,cap,"winding$(charge)_penalty_pass$(pass)"))
        previousS=S
        passed && return result
    end
    error("Sector convergence not achieved; all stages and checkpoints preserved in $stem.json")
end
if abspath(PROGRAM_FILE)==@__FILE__
    @assert length(ARGS) in [5,6,7] "family L w N_up q [maxcap] [resume-checkpoint]"
    run_sector_penalty(ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3]),parse(Int,ARGS[4]),parse(Int,ARGS[5]);
        maxcap=length(ARGS)>=6 ? parse(Int,ARGS[6]) : 256,resume=length(ARGS)==7 ? ARGS[7] : nothing)
end
