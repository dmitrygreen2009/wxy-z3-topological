isdefined(Main,:cylinder) || include("../src/model.jl")
using ITensorInfiniteMPS, KrylovKit, Serialization
include(joinpath(pkgdir(ITensorInfiniteMPS),"examples","vumps","src","vumps_subspace_expansion.jl"))

include("../src/infinite_model.jl")
include("../src/vumps_audit.jl")
include("../src/transfer_analysis.jl")

function measure_infinite(psi,H,family,w,cap,iteration;tag="",ordering="matter_first")
    measurement_started=time()
    n=length(psi.AL);@assert n%(9w)==0
    cell_slices=n÷(9w)
    slice_cuts=collect(9w:9w:n)
    stem="results/infinite_$(family)_w$(w)_chi$(cap)"*tag
    energies=expect(psi,H)
    # Densify only for measuring all transfer charge sectors, not optimizing.
    # A zero-flux block alone can miss the leading non-neutral eigenvalue.
    is_qn=hasqns(siteind(psi.AL,1))
    A=is_qn ? InfiniteMPS([dense(psi.AL[j]) for j=1:n],translator(psi.AL)) : psi.AL
    T=TransferMatrix(A)
    spectrum=transfer_spectrum(A)
    vals=spectrum["values"];xi=spectrum["xi_cells"]
    ratios=abs.(vals./vals[1])
    entropies=[]
    for j=1:n
        p=abs2.(svdvals(Array(psi.C[j],inds(psi.C[j])...)));p./=sum(p)
        push!(entropies,-sum(x>0 ? x*log(x) : 0.0 for x in p))
    end
    raw_errors=[];phase_errors=[];left_isometry_errors=[];right_isometry_errors=[];center_norm_errors=[]
    for j=1:n
        A=psi.AL[j]*psi.C[j];B=psi.C[j-1]*psi.AR[j]
        overlap=inner(A,B)
        push!(raw_errors,norm(A-B))
        push!(phase_errors,norm(A-(abs(overlap)>0 ? conj(overlap)/abs(overlap) : 1)*B))
        r=commonind(psi.AL[j],psi.AL[j+1]);l=commonind(psi.AR[j-1],psi.AR[j])
        push!(left_isometry_errors,norm(dense(psi.AL[j]*dag(prime(psi.AL[j],r)))-dense(delta(r,dag(prime(r)))))/sqrt(dim(r)))
        push!(right_isometry_errors,norm(dense(psi.AR[j]*dag(prime(psi.AR[j],l)))-dense(delta(l,dag(prime(l)))))/sqrt(dim(l)))
        push!(center_norm_errors,abs(norm(psi.C[j])-1))
    end
    sz_profile=[real(expect(psi,"Sz",j)) for j=1:n]
    r=Dict("family"=>family,"width"=>w,"chi"=>maxlinkdim(psi),"cap"=>cap,"iteration"=>iteration,"cell_spins"=>n,
        "energy_cell"=>real(sum(energies)),"entropy_bonds"=>entropies,"spatial_entropy"=>entropies[end],
        "spatial_entropies_at_slice_cuts"=>entropies[slice_cuts],
        "slice_cut_entropy_modulation"=>maximum(entropies[slice_cuts])-minimum(entropies[slice_cuts]),
        "spatial_entropy_slice_mean"=>sum(entropies[slice_cuts])/cell_slices,
        "mean_sz"=>sum(sz_profile)/n,"site_sz_profile_infinite_order"=>sz_profile,"ordering_version"=>2,
        "mean_abs_splus"=>is_qn ? 0.0 : sum(abs(expect(psi,"S+",j)) for j=1:n)/n,
        "u1_conserving_ansatz"=>is_qn,"measurement_tag"=>tag,"infinite_ordering"=>ordering,
        "transfer_converged_eigenpairs"=>spectrum["converged"],"transfer_eigenvalue_magnitudes"=>abs.(vals),
        "transfer_eigenvalues_real"=>real.(vals),"transfer_eigenvalues_imag"=>imag.(vals),
        "xi_cells"=>isfinite(xi) ? xi : "Inf","xi_slices"=>isfinite(xi) ? cell_slices*xi : "Inf",
        "xi_physical_axial"=>isfinite(xi) ? cell_slices*xi*(family=="zigzag" ? 1.5 : sqrt(3)/2) : "Inf",
        "correlation_length_units"=>"Cell size is recorded in cell_slices; physical axial lengths use nearest-neighbor honeycomb distance one",
        "transfer_residuals"=>spectrum["residuals"],
        "fixed_point_rank_detected"=>spectrum["fixed_point_rank_detected"],
        "fixed_point_gram_eigenvalues"=>spectrum["fixed_point_gram_eigenvalues"],
        "fixed_point_residuals"=>spectrum["fixed_point_residuals"],
        "transfer_gap_above_residual_scale"=>spectrum["ratio_gap_above_residual_scale"],
        "transfer_normalization_error"=>abs(abs(vals[1])-1),
        "subleading_transfer_above_noise"=>length(vals)>=2 && ratios[2]>1e-12,
        "canonical_error_raw"=>maximum(raw_errors),"canonical_error"=>maximum(phase_errors),
        "left_isometry_error"=>maximum(left_isometry_errors),"right_isometry_error"=>maximum(right_isometry_errors),
        "center_normalization_error"=>maximum(center_norm_errors),
        "process_peak_rss_bytes"=>Sys.maxrss())
    # Species-one Bell dimers use leg one on A and leg two on B,
    # all distinct physical gauges. Each contributes -1/sqrt(3), all
    # other exchanges average to zero, and half filling is possible.
    r["disjoint_dimer_variational_upper_bound_cell"]=-2cell_slices*w/sqrt(3)
    # Alternatively every A star and its three physical gauges form a
    # disjoint validated six-spin cluster. Fixed-number B matter gives
    # zero mean to the remaining exchanges. Round the trial bound upward.
    r["isolated_A_star_variational_upper_bound_cell"]=-2.4030921cell_slices*w
    r["energy_above_known_trial_state"]=r["energy_cell"]>-2.4030921cell_slices*w+1e-9
    r["variational_bound_interpretation"]="Energy above the explicit trial state excludes a ground-state candidate irrespective of projected solver residual; satisfying this bound does not certify convergence."
    r["kind"]="infinite"
    r["validation_fixture"]=occursin("_measurement_validation",tag)
    r["translation_period_restriction"]="$(cell_slices) axial slices"
    r["topological_flux_or_MES_identified"]=false
    r["number_background"]=infinite_number_background(psi)
    r["bulk_ground_filling_certified"]=false
    r["physical_bulk_gap_certified"]=false
    r["fixed_point_detection_scope"]="Detected multiplicity is a lower bound from two independent starts, not a complete peripheral-spectrum certificate."
    r["physical_circumference"]=family=="zigzag" ? sqrt(3)*w : 3.0w
    r["cell_slices"]=cell_slices;r["physical_spins"]=n;r["energy_per_vertex"]=r["energy_cell"]/(2cell_slices*w)
    r["measurement_runtime_seconds"]=time()-measurement_started
    audit_path=stem*"_solver_audit.json"
    if isfile(audit_path)
        a=JSON3.read(read(audit_path,String),Dict{String,Any})
        r["audit"]=a;r["git_commit"]=a["git_commit"]
        r["number_background"]=get(a,"number_background",r["number_background"])
        r["measurement_tensors_have_qns"]=is_qn
        r["u1_conserving_ansatz"]=get(a,"optimization_u1_conserving_ansatz",is_qn)
        r["runtime_seconds"]=a["runtime_seconds"]+r["measurement_runtime_seconds"]
        r["solver_residual"]=a["canonical_solver_residual"]
        if r["u1_conserving_ansatz"] && r["number_background"]["label"]=="unrestricted"
            r["number_background"]=Dict("label"=>"U1_background_unrecorded","density_enforced"=>nothing,"interpretation"=>"Original U1 ansatz was recorded, but legacy dense measurement tensors do not retain the physical QN background. Do not infer half filling from zero flux.")
        end
        if haskey(a,"source_spatial_entropy")
            r["entropy_change_recanonicalization"]=r["spatial_entropy"]-a["source_spatial_entropy"]
            r["energy_change_recanonicalization"]=r["energy_cell"]-a["source_energy_cell"]
        end
    else
        r["git_commit"]=nothing
        r["historical_audit_gap"]="Original run predates launch provenance; measurement and legacy logs are preserved."
    end
    target_density=get(r["number_background"],"physical_number_density",nothing)
    r["number_density_error_vs_background"]=target_density===nothing ? nothing : abs(r["mean_sz"]+.5-target_density)
    sector=r["number_background"]["label"]
    Ly=round(r["physical_circumference"];digits=5)
    measurement_seed=get(get(r,"audit",Dict()),"random_seed",7103)
    final_path="results/checkpoints/infinite_$(family)_w$(w)_Ly$(Ly)_cell$(cell_slices)_N$(n)_chi$(maxlinkdim(psi))_$(sector)_seed$(measurement_seed)_v$(RUN_FORMAT_VERSION)$(tag)_stage$(iteration)_complete.jls"
    manifest=save_checkpoint(final_path,psi,r)
    atomic_json(stem*".json",manifest)
    open("results/infinite_$(family)_w$(w)$(tag)_history.jsonl","a") do io;println(io,JSON3.write(r));end
    println("MEASURED: ",JSON3.write(r));flush(stdout)
end

function run_infinite(family,w,chi;ordering="matter_first",tag="")
    Random.seed!(7103);BLAS.set_num_threads(1)
    n=18w
    initstate(j)=isodd(j) ? "Up" : "Dn"
    s=infsiteinds("S=1/2",n;conserve_qns=false,initstate)
    psi=InfMPS(s,initstate)
    for j=1:n
        theta=0.2+rand()*1.1;phase=cis(2pi*rand())
        rotation=ComplexF64[cos(theta) -conj(phase)*sin(theta);phase*sin(theta) cos(theta)]
        u=ITensor(rotation,s[j]',dag(s[j]))
        psi.AL[j]=noprime(u*psi.AL[j]);psi.AR[j]=noprime(u*psi.AR[j])
    end
    H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering),s)
    caps=unique(vcat([min(c,chi) for c in [2,4,8,16,32,64,128]],chi))
    for (iteration,cap) in enumerate(caps)
        estimate_memory(psi,cap;label="$(family) infinite w$(w) expansion")
        psi=subspace_expansion(psi,H;cutoff=1e-9,maxdim=cap)
        psi=audited_vumps(H,psi;family,w,cap,ordering,tag,tol=1e-7,maxiter=30,
            solver_tol=x->max(x/10,1e-10),solver_tolerance_rule="max(residual/10,1e-10)")
        measure_infinite(psi,H,family,w,cap,iteration;tag,ordering)
    end
    psi=audited_vumps(H,psi;family,w,cap=chi,ordering,tag,tol=1e-7,maxiter=40,
        solver_tol=x->max(x/10,1e-10),solver_tolerance_rule="max(residual/10,1e-10)")
    measure_infinite(psi,H,family,w,chi,length(caps)+1;tag,ordering)
end
if abspath(PROGRAM_FILE)==@__FILE__
    run_infinite(ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3]))
end
