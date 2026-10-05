# Exact microscopic, individual-site infinite optimization. No effective theory.
include("infinite.jl")
include("../src/infinite_winding_basis.jl")
include("../src/infinite_block_initializer.jl")
include("../src/canonical_analysis.jl")
include("../src/transfer_uniqueness.jl")
const INFINITE_WINDING_DRIVER_SHA=bytes2hex(sha256(read(@__FILE__)))

function validate_infinite_winding_state(psi,layout,nup)
    support=infinite_winding_boundary_support(psi,layout);sites=siteinds(only,psi);loops=[]
    for group in eachindex(layout.charges),offset in [0,layout.n]
        z=expect(psi,cgs_product_mpo(sites,infinite_winding_diagonal_spec(layout,group;offset)))
        error=abs(z-cis(2pi*layout.charges[group]/3))
        @assert error<1e-10 && abs(1-abs2(z))<1e-10
        push!(loops,Dict("group"=>group,"offset_sites"=>offset,"charge"=>layout.charges[group],
            "expectation_real"=>real(z),"expectation_imag"=>imag(z),"variance"=>1-abs2(z),"purity_error"=>error))
    end
    background=infinite_number_background(psi)
    @assert background["density_enforced"] && background["target_average_nup_per_cell"]==nup
    Dict("boundary_support"=>support,"individual_winding_loops"=>loops,"number_background"=>background)
end

function run_infinite_winding_qn(family,w,nup,charge,maxcap;cell_slices=2,seed=7260,resume=nothing,maxiter=60,validation_fixture=false)
    started=time();step=family=="zigzag" ? 1 : 2
    @assert cell_slices%step==0
    charges=fill(charge,cell_slices÷step);layout=infinite_winding_cell_layout(family,w,cell_slices,charges)
    geometry="geometry/infinite_$(family)_w$(w)$(cell_slices==2 ? "" : "_cell$(cell_slices)")_star.json"
    @assert isfile(geometry) && JSON3.read(read(geometry,String))["cell_spins"]==layout.n
    if resume===nothing
        phi,layout=infinite_winding_initial_cell(family,w,nup,charges;cell_slices,seed,linkdims=2)
        psi=infinite_cell_product(phi);initial_cap=2;initialization="Exact finite-cell product with recorded number and individual winding QNs"
    else
        actual=resolve_checkpoint(resume);meta=JSON3.read(read(replace(actual,".jls"=>".json"),String),Dict{String,Any})
        @assert checkpoint_basis(meta)=="exact_matter_charge_basis"
        @assert meta["family"]==family && meta["width"]==w && meta["cell_slices"]==cell_slices
        settings=get(meta,"solver_settings",get(get(meta,"audit",Dict()),"solver_settings",Dict()))
        @assert Int.(settings["individual_loop_charges"])==charges
        @assert [Int.(row) for row in settings["onsite_winding_weights"]]==[collect(layout.weights[j,:]) for j=1:layout.n]
        psi=load_state(actual);initial_cap=Int(meta["cap"]);initialization=actual
        @assert hasqns(siteind(psi.AL,1)) "Resume the QN optimization checkpoint, not a dense measurement copy"
    end
    @assert maximum_bond_dimension(psi)<=max(initial_cap,2)<=maxcap
    validate(state)=validate_infinite_winding_state(state,layout,nup)
    validate(psi)
    caps=sort(unique(vcat(initial_cap,[c for c in [2,4,8,16,32,64,128,256] if initial_cap<c<maxcap],maxcap)))
    tag="_wind$(join(charges,"x"))_Nup$(nup)_charge_basis_star_seed$(seed)_run$(Dates.format(now(UTC),"yyyymmddTHHMMSSsss"))"
    previous=nothing;stage_results=[]
    for (stage,cap) in enumerate(caps)
        estimate_memory(psi,cap;label="Exact individual-winding infinite $family w$w density $(nup/layout.n)")
        H=InfiniteSum{MPO}(infinite_real_charge_opsum(family,w;cell_slices),siteinds(only,psi))
        psi=subspace_expansion(psi,H;cutoff=1e-12,maxdim=cap);validate(psi)
        audit=run_provenance(;seed,solver="Official ITensorInfiniteMPS with individual microscopic winding-boundary blocks",
            settings=Dict("basis"=>"exact_matter_charge_basis","hamiltonian_representation"=>"real_closed_form",
                "individual_loop_charges"=>charges,"winding_field_names"=>layout.field_names,
                "onsite_winding_weights"=>[collect(layout.weights[j,:]) for j=1:layout.n],
                "winding_background_offsets"=>[collect(layout.background[j,:]) for j=1:layout.n],
                "winding_cell_boundary_bonds"=>layout.ends,"physical_nup_per_cell"=>nup,
                "solver_tolerance"=>1e-7,"local_eigensolver_tolerance"=>1e-12,"maxiter"=>maxiter,
                "subspace_expansion_cutoff"=>1e-12,"loop_tolerance"=>1e-10,"update_algorithm"=>"sequential"),
            initialization,conserved_quantum_numbers=vcat(["Fixed U1 mean N_up=$nup per $(layout.n) sites"],["Individual loop $(g) charge=$(q)" for (g,q) in enumerate(charges)]))
        audit["executed_driver_sha256"]=INFINITE_WINDING_DRIVER_SHA
        audit["optimization_u1_conserving_ansatz"]=true
        psi=audited_vumps(H,psi;family,w,cap,maxiter,tol=1e-7,solver_tol=x->1e-12,
            ordering="star",tag,seed,audit,update_algorithm="sequential",state_validator=validate)
        validation=validate(psi)
        number_background=infinite_number_background(psi);sector=number_background["label"]
        latest="results/checkpoints/infinite_$(family)_w$(w)_cell$(cell_slices)_N$(layout.n)_chi$(cap)_$(sector)_seed$(seed)_v$(RUN_FORMAT_VERSION)$(tag)_latest.jls"
        source_sha=open(io->bytes2hex(sha256(io)),latest)
        # Correct Schmidt centers on the same AL state only for a unique fixed point.
        raw_error=maximum(begin
            a=psi.AL[j]*psi.C[j];b=psi.C[j-1]*psi.AR[j];z=inner(a,b)
            norm(a-(abs(z)>0 ? conj(z)/abs(z) : 1)*b)
        end for j=1:layout.n)
        measurement=psi;canonical_action="Original consistent canonical state";uniqueness=nothing
        if raw_error>1e-10
            A=InfiniteMPS([dense(psi.AL[j]) for j=1:layout.n],translator(psi.AL))
            uniqueness=transfer_uniqueness_certificate(A)
            if uniqueness["primitive_peripheral_spectrum_certified"]
                measurement,_=canonicalize_left(psi.AL)
                canonical_action="Same AL state recanonicalized for measurement; optimization checkpoint unchanged"
            else
                canonical_action="Primitive fixed point not certified: original centers retained and entropy consistency unresolved"
            end
        end
        MH=InfiniteSum{MPO}(infinite_real_charge_opsum(family,w;cell_slices),siteinds(only,measurement))
        captured=Ref{Any}()
        function interpret(r)
            profile=r["site_sz_profile_infinite_order"]
            physical_profile=copy(profile)
            for row=0:cell_slices-1,j=0:w-1,first in [9w*row+9j+1,9w*row+9j+6]
                physical_profile[first:first+2].=sum(profile[first:first+2])/3
            end
            r["mode_sz_profile_exact_charge_basis"]=profile
            r["site_sz_profile_infinite_order"]=physical_profile
            r["physical_matter_sz_derivation"]="A definite exact loop state is invariant in expectation under its nontrivial cyclic matter permutation. Each species has one third of basis-invariant star magnetization."
            r["physical_density_profile_consistent"]=r["canonical_error"]<1e-10
            r["basis"]="exact_matter_charge_basis";r["physical_geometry"]=geometry
            r["validation_fixture"]=validation_fixture
            r["individual_winding_validation"]=validation;r["individual_loop_charges"]=charges
            r["unit_cell_nup"]=nup;r["unit_cell_sz"]=nup-layout.n/2
            r["measured_unit_cell_nup"]=layout.n/2+sum(physical_profile)
            r["filling_fraction"]=nup/layout.n;r["filling_selection"]="Fixed rational U1 background, not thermodynamic variational selection"
            r["measured_number_consistent_with_background"]=abs(r["measured_unit_cell_nup"]-nup)<1e-9
            if r["canonical_error"]<1e-10
                @assert r["measured_number_consistent_with_background"]
            end
            r["measurement_canonical_action"]=canonical_action;r["optimization_canonical_error_before_measurement"]=raw_error
            r["measurement_transfer_uniqueness_certificate"]=uniqueness
            r["optimization_resume_checkpoint"]=latest;r["optimization_checkpoint_sha256"]=source_sha
            r["individual_microscopic_loop_sectors_verified"]=true
            r["topological_flux_or_MES_identified"]=false;r["eligible_for_topological_entropy_fit"]=false
            r["within_winding_pattern_exclusion_supported_by_trial"]=false
            r["fixed_density_exclusion_scope"]="Comparison across all loop sectors at this mean number; not a within-pattern bound"
            r["entropy_measurement_consistent"]=r["canonical_error"]<1e-10
            r["runtime_seconds_total_driver"]=time()-started
            captured[]=r
        end
        measure_infinite(measurement,MH,family,w,cap,stage;tag,ordering="star",result_transform=interpret)
        @assert open(io->bytes2hex(sha256(io)),latest)==source_sha
        r=captured[]
        quality=all(r[k]<1e-10 for k in ["canonical_error","left_isometry_error","right_isometry_error","center_normalization_error","transfer_normalization_error"])
        differences=previous===nothing ? nothing : Dict("energy_per_vertex"=>abs(r["energy_per_vertex"]-previous["energy_per_vertex"]),
            "spatial_entropy"=>abs(r["spatial_entropy"]-previous["spatial_entropy"]))
        push!(stage_results,Dict("cap"=>cap,"actual_bond_dimension"=>r["chi"],"energy_cell"=>r["energy_cell"],
            "entropy"=>r["spatial_entropy"],"xi_cells"=>r["xi_cells"],"solver_residual"=>r["solver_residual"],
            "measurement_consistent"=>quality,"change_previous_cap"=>differences,"resume_checkpoint"=>latest))
        atomic_json("results/infinite_$(family)_w$(w)$(tag)_bond_study.json",Dict("records"=>stage_results,
            "nup_per_cell"=>nup,"sz_per_cell"=>nup-layout.n/2,"filling"=>nup/layout.n,"charges"=>charges,
            "interpretation"=>"Fixed-number-background and disjoint-loop-pattern study; MES, thermodynamic filling and 2D phase remain uncertified."))
        previous=r
    end
    stage_results
end
if abspath(PROGRAM_FILE)==@__FILE__
    @assert length(ARGS) in [5,6]
    run_infinite_winding_qn(ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3]),parse(Int,ARGS[4]),parse(Int,ARGS[5]);resume=length(ARGS)==6 ? ARGS[6] : nothing)
end
