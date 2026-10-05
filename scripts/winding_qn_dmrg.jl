# Exact microscopic Hamiltonian in a number-preserving local matter basis.
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/matter_charge_basis.jl")
include("../src/cgs_operator_mpo.jl")
const WINDING_QN_DRIVER_SHA=bytes2hex(sha256(read(@__FILE__)))
function run_winding_qn(family,L,w,nup,charge;seed=7254,maxcap=512,resume=nothing,krylovdim=12,eigsolve_maxiter=30,noise=[1e-5,1e-6,1e-7,0.0])
    Random.seed!(seed);started=time();ordering="star";lat=cylinder(family,L,w;ordering)
    table="geometry/cgs_cycles/$(family)_L$(L)_w$(w)_star.json"
    defs=JSON3.read(read(table,String),Dict{String,Any});cycle=first(c for c in defs["cycles"] if c["winding_number"]==1)
    g=JSON3.read(read(defs["parent_geometry"],String),Dict{String,Any})
    @assert g["physical_spins"]==lat.n && Int.(g["mps_order_physical_spin_ids"])==lat.order && g["spatial_cut_mps_bond"]==lat.cut
    p=Int.(cycle["edge_exponents"]);weights=winding_site_weights(lat,p)
    if resume===nothing
        sites=winding_qn_sites(weights);state=winding_initial_state(weights,nup,charge)
        psi=random_mps(ComplexF64,sites,state;linkdims=4)
    else
        psi=load_state(resume);sites=siteinds(psi)
        @assert length(psi)==lat.n
        previous=JSON3.read(read(replace(resume,".jls"=>".json"),String),Dict{String,Any})
        settings=get(previous,"solver_settings",get(get(previous,"audit",Dict()),"solver_settings",Dict()))
        @assert settings["basis"]=="exact_matter_charge_basis" && Int.(settings["onsite_winding_weights"])==weights
    end
    @assert val(flux(psi),"Sz")==2nup-lat.n && mod(val(flux(psi),"Winding"),3)==charge
    H=MPO(charge_basis_opsum(lat),sites)
    U=cgs_product_mpo(sites,(;triplets=Tuple{Int,Int,Int}[],gauges=collect(enumerate(weights)),first_site=1,last_site=lat.n))
    known_energy=nothing
    for edpath in ["results/sector_penalty_ed_audit.json","results/winding_physical_sector_ed.json"]
        isfile(edpath) || continue
        ed=JSON3.read(read(edpath,String),Dict{String,Any})
        for row in ed["records"]
            if row["family"]==family && row["L"]==L && row["width"]==w && row["nup"]==nup
                matches=[sector["energy"] for sector in row["sectors"] if sector["charge"]==charge]
                isempty(matches) || (known_energy=only(matches))
            end
        end
    end
    audit=run_provenance(;seed,solver="ITensorMPS exact U1 x Z3 winding-QN DMRG",
        settings=Dict("basis"=>"exact_matter_charge_basis","onsite_winding_weights"=>weights,"winding_charge"=>charge,
            "edge_exponents"=>p,"matter_hilbert_space_dimension"=>8,"local_roundoff_cutoff"=>1e-14,
            "cutoff"=>1e-13,"sweeps_per_stage"=>12,"krylovdim"=>krylovdim,"eigsolve_maxiter"=>eigsolve_maxiter,"noise_schedule"=>noise,"noise_application"=>"First stage only; subsequent stages have zero noise to avoid rotating degenerate ground states","eigsolve_tol"=>1e-13,
            "purity_tolerance"=>1e-10,"variance_tolerance"=>1e-8,"energy_drift_tolerance"=>1e-9,"entropy_drift_tolerance"=>1e-5,"truncation_error_window"=>"Last four zero-noise sweeps",
            "known_fixed_sector_ed_energy"=>known_energy,"known_energy_tolerance"=>1e-9),
        initialization=resume===nothing ? "Random MPS in exact number/winding QN block" : resume,
        conserved_quantum_numbers=["Physical N_up=$nup","Exact microscopic winding charge=$charge"])
    audit["executed_driver_sha256"]=WINDING_QN_DRIVER_SHA
    stem="results/$(family)_L$(L)_w$(w)_Nup$(nup)_winding$(charge)_qn_seed$(seed)"
    records=[];previousS=nothing
    for cap in unique(min.(maxcap,[32,64,128,256,512,1024])),pass in 1:2
        estimate_memory(psi,cap;label="Direct winding-QN $family L$L w$w q$charge")
        stage_noise=isempty(records) ? noise : [0.0]
        obs=finite_checkpoint_observer(lat;family,L,w,cap,seed,stage="direct_winding$(charge)_pass$(pass)",ordering,audit,cutoff=1e-13,noise=stage_noise,nup)
        _,psi=dmrg(H,psi;nsweeps=12,maxdim=cap,cutoff=1e-13,noise=stage_noise,eigsolve_krylovdim=krylovdim,eigsolve_maxiter,eigsolve_tol=1e-13,observer=obs,outputlevel=1)
        @assert val(flux(psi),"Sz")==2nup-lat.n && mod(val(flux(psi),"Winding"),3)==charge
        E=real(inner(psi',H,psi)/inner(psi,psi));variance=real(inner(H,psi,H,psi)/inner(psi,psi))-E^2
        @assert variance>=-1e-8
        mean=inner(psi',U,psi)/inner(psi,psi);purity=abs(mean-cis(2pi*charge/3))
        S,probs=entropy_at(psi,lat.cut);es=energies(obs);drift=abs(es[end]-es[end-1]);entropy_drift=previousS===nothing ? nothing : abs(S-previousS)
        passed=purity<1e-10 && variance<1e-8 && drift<1e-9 && maximum(truncerrors(obs)[end-3:end])<1e-10 && entropy_drift!==nothing && entropy_drift<1e-5 && (known_energy===nothing || abs(E-known_energy)<1e-9)
        push!(records,Dict("cap"=>cap,"pass"=>pass,"stage_noise_schedule"=>stage_noise,"energy"=>E,"entropy"=>S,"maxlinkdim"=>maxlinkdim(psi),
            "schmidt_probabilities"=>probs,"sweep_energies"=>es,"sweep_max_truncation_errors"=>truncerrors(obs),
            "energy_variance"=>variance,"loop_expectation_real"=>real(mean),"loop_expectation_imag"=>imag(mean),
            "loop_variance"=>1-abs2(mean),"sector_purity_error"=>purity,"entropy_change_previous_stage"=>entropy_drift,
            "known_fixed_sector_ed_energy"=>known_energy,"energy_error_vs_sector_ed"=>known_energy===nothing ? nothing : E-known_energy,
            "converged_within_fixed_number_and_loop_sector"=>passed))
        result=Dict("family"=>family,"length"=>L,"width"=>w,"physical_circumference"=>lat.circumference,"spins"=>lat.n,"vertices"=>lat.nv,
            "ordering"=>ordering,"seed"=>seed,"basis"=>"exact_matter_charge_basis","nup"=>nup,"total_sz"=>nup-lat.n/2,
            "filling_fraction"=>nup/lat.n,"filling_selection"=>"Fixed physical number and winding MPS QNs",
            "cut"=>lat.cut,"loop_label"=>cycle["label"],"geometry_cycle_definitions"=>table,"physical_geometry"=>defs["parent_geometry"],
            "records"=>records,"audit"=>audit,"git_commit"=>audit["git_commit"],"runtime_seconds"=>time()-started,"correlation_length"=>nothing,
            "entropy_interpretation"=>"Full-state entropy invariant under local matter basis changes wholly inside each spatial partition",
            "interpretation"=>"Exact microscopic fixed-number/winding sector; no Hilbert-space truncation beyond MPS bond truncation, no MES or 2D phase certificate.")
        atomic_json(stem*".json",completed_finite_checkpoint(psi,result,cap,"direct_winding$(charge)_pass$(pass)"))
        previousS=S;passed && return result
    end
    error("Direct sector-QN convergence incomplete; saved data and checkpoints preserved")
end
if abspath(PROGRAM_FILE)==@__FILE__
    @assert length(ARGS) in [5,6,7]
    run_winding_qn(ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3]),parse(Int,ARGS[4]),parse(Int,ARGS[5]);maxcap=length(ARGS)>=6 ? parse(Int,ARGS[6]) : 512,resume=length(ARGS)==7 ? ARGS[7] : nothing)
end
