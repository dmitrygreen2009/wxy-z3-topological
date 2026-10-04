isdefined(Main,:cylinder) || include("../src/model.jl")
using Serialization
Random.seed!(7103); BLAS.set_num_threads(1)
function run_cylinder(family,L,w,chi;seed=7103,ordering="axial",nup=nothing)
Random.seed!(seed)
started=time()
lat=cylinder(family,L,w;ordering)
sector_nup=nup===nothing ? cld(lat.n,2) : nup
@assert 0<=sector_nup<=lat.n
audit=run_provenance(;seed,solver="ITensorMPS finite two-site DMRG",
    settings=Dict("target_chi"=>chi,"stage_sweeps"=>8,"refinement_sweeps"=>6,"cutoff"=>1e-10,"refinement_cutoff"=>1e-11,
        "sweep_energy_absolute_tolerance"=>1e-8,"bond_energy_per_vertex_tolerance"=>1e-6,
        "bond_entropy_tolerance"=>1e-3,"refinement_entropy_tolerance"=>1e-4,"truncation_error_target"=>1e-8,"length_entropy_tolerance"=>1e-3,
        "krylovdim"=>8,"refinement_krylovdim"=>12,"convergence_criteria"=>"Recorded energy, entropy, truncation, chi and length drift; sweep completion alone is not convergence"),
    initialization="Random QN MPS at N_up=$sector_nup from an adjusted alternating pattern, initial link dimension 8",conserved_quantum_numbers=["total N_up=$sector_nup"])
mkpath("geometry");export_cylinder_geometry(family,L,w;ordering)
stem="results/$(family)_L$(L)_w$(w)_chi$(chi)"*(ordering=="axial" ? "" : "_"*ordering)*(seed==7103 ? "" : "_seed$(seed)")
nup===nothing || (stem*="_Nup$(sector_nup)")
sites=siteinds("S=1/2",lat.n;conserve_qns=true)
H=MPO(lat.os,sites)
state=[isodd(i) ? "Up" : "Dn" for i=1:lat.n]
difference=sector_nup-count(==("Up"),state)
if difference!=0
    source=difference>0 ? "Dn" : "Up";destination=difference>0 ? "Up" : "Dn"
    candidates=shuffle(findall(==(source),state))
    state[candidates[1:abs(difference)]].=destination
end
@assert count(==("Up"),state)==sector_nup
psi=random_mps(sites,state;linkdims=8)
@assert val(flux(psi),"Sz")==2sector_nup-lat.n
records=[]
for cap in unique([min(32,chi),min(64,chi),min(128,chi),min(256,chi),chi])
    estimate_memory(psi,cap;label="$(family) L$(L) w$(w) DMRG")
    obs=finite_checkpoint_observer(lat;family,L,w,cap,seed,stage="cap$(cap)",ordering,audit,nup=sector_nup,cutoff=1e-10,noise=[1e-5,1e-6,1e-7,0])
    e,psi=dmrg(H,psi;nsweeps=8,maxdim=cap,cutoff=1e-10,
        noise=[1e-5,1e-6,1e-7,0],eigsolve_krylovdim=8,outputlevel=1,observer=obs)
    S,p=entropy_at(psi,lat.cut)
    push!(records,Dict("cap"=>cap,"energy"=>e,"entropy"=>S,"maxlinkdim"=>maxlinkdim(psi),"schmidt_probabilities"=>p,
        "sweep_energies"=>energies(obs),"sweep_max_truncation_errors"=>truncerrors(obs),"process_peak_rss_bytes"=>Sys.maxrss()))
    result=Dict("family"=>family,"length"=>L,"width"=>w,"physical_circumference"=>lat.circumference,
        "spins"=>lat.n,"vertices"=>lat.nv,"nup"=>sector_nup,"explicit_number_sector"=>nup!==nothing,"ground_state_scope"=>"Fixed number sector; global ground sector requires an independent search","cut"=>lat.cut,"records"=>records,"seed"=>seed,"ordering"=>ordering,
        "legs"=>lat.legs,"vertices_table"=>lat.verts,"mps_order"=>lat.order,
        "audit"=>audit,"git_commit"=>audit["git_commit"],"runtime_seconds"=>time()-started,
        "correlation_length"=>nothing,"correlation_length_status"=>"Use separately optimized infinite-state transfer spectrum")
    atomic_json(stem*".json",result)
    manifest=completed_finite_checkpoint(psi,result,cap,"stage")
    atomic_json(stem*".json",manifest)
end
Sbefore,_=entropy_at(psi,lat.cut)
obs=finite_checkpoint_observer(lat;family,L,w,cap=chi,seed,stage="refinement",ordering,audit,nup=sector_nup,cutoff=1e-11,noise=0)
e,psi=dmrg(H,psi;nsweeps=6,maxdim=chi,cutoff=1e-11,noise=0,
    eigsolve_krylovdim=12,outputlevel=1,observer=obs)
S,p=entropy_at(psi,lat.cut)
push!(records,Dict("cap"=>chi,"energy"=>e,"entropy"=>S,"maxlinkdim"=>maxlinkdim(psi),
    "schmidt_probabilities"=>p,"sweep_energies"=>energies(obs),"sweep_max_truncation_errors"=>truncerrors(obs),
    "entropy_change_refinement"=>S-Sbefore,"process_peak_rss_bytes"=>Sys.maxrss()))
result=Dict("family"=>family,"length"=>L,"width"=>w,"physical_circumference"=>lat.circumference,
    "spins"=>lat.n,"vertices"=>lat.nv,"nup"=>sector_nup,"explicit_number_sector"=>nup!==nothing,"ground_state_scope"=>"Fixed number sector; global ground sector requires an independent search","cut"=>lat.cut,"records"=>records,"seed"=>seed,"ordering"=>ordering,
    "legs"=>lat.legs,"vertices_table"=>lat.verts,"mps_order"=>lat.order,
        "audit"=>audit,"git_commit"=>audit["git_commit"],"runtime_seconds"=>time()-started,
        "correlation_length"=>nothing,"correlation_length_status"=>"Use separately optimized infinite-state transfer spectrum")
atomic_json(stem*".json",result)
manifest=completed_finite_checkpoint(psi,result,chi,"refined")
atomic_json(stem*".json",manifest)
end
if abspath(PROGRAM_FILE)==@__FILE__
    run_cylinder(ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3]),parse(Int,ARGS[4]);ordering=length(ARGS)>=5 ? ARGS[5] : "axial")
end
