include("cylinders.jl")
include("../src/cgs.jl")
function run_projected_zigzag(L,chi;seed=7112)
    Random.seed!(seed);started=time();family="zigzag";w=1;ordering="star"
    lat=cylinder(family,L,w;ordering);sites=siteinds("S=1/2",lat.n;conserve_qns=true)
    H=MPO(lat.os,sites);psi=random_mps(ComplexF64,sites,[isodd(i) ? "Up" : "Dn" for i=1:lat.n];linkdims=8)
    audit=run_provenance(;seed,solver="ITensorMPS DMRG with exact microscopic CGS-sector initialization",
        settings=Dict("chi"=>chi,"cutoff"=>1e-12,"projection_cutoff"=>1e-13,"projection_maxdim"=>256,"sweeps_per_cap"=>10,"noise"=>0,"krylovdim"=>12,"selected_local_loop_charge"=>0),
        initialization="Complex random QN MPS followed by exact P0=(I+U+U^2)/3 for every parallel-edge two-cycle",conserved_quantum_numbers=["total Sz","initial trivial parallel-edge cycle charges"])
    projection_records=[]
    for t=0:L-1
        v=findfirst(==((t,0,0)),lat.verts);p=zeros(Int,lat.n-3lat.nv)
        p[lat.legs[v][1]]=1;p[lat.legs[v][3]]=2
        psi,record=cgs_project(psi,lat,p;maxdim=256,cutoff=1e-13)
        record["slice"]=t;record["edge_exponents"]=p
        @assert record["purity_error"]<1e-8 "CGS projection purity failed; do not relax tolerance"
        push!(projection_records,record)
    end
    stem="results/zigzag_L$(L)_w1_chi$(chi)_star_projected0_seed$(seed)"
    records=[]
    for cap in unique(vcat([min(c,chi) for c in [32,64,128,256,512]],chi))
        estimate_memory(psi,cap;label="Projected zigzag L$(L) width1 DMRG")
        obs=finite_checkpoint_observer(lat;family,L,w,cap,seed,stage="projected0_cap$(cap)",ordering,audit,cutoff=1e-12,noise=0)
        E,psi=dmrg(H,psi;nsweeps=10,maxdim=cap,cutoff=1e-12,noise=0,eigsolve_krylovdim=12,observer=obs,outputlevel=1)
        S,p=entropy_at(psi,lat.cut)
        push!(records,Dict("cap"=>cap,"energy"=>E,"entropy"=>S,"schmidt_probabilities"=>p,"maxlinkdim"=>maxlinkdim(psi),"sweep_energies"=>energies(obs),"sweep_max_truncation_errors"=>truncerrors(obs)))
        result=Dict("family"=>family,"width"=>w,"length"=>L,"physical_circumference"=>lat.circumference,"spins"=>lat.n,"vertices"=>lat.nv,
            "ordering"=>ordering,"seed"=>seed,"cut"=>lat.cut,"records"=>records,"projection_records"=>projection_records,"audit"=>audit,
            "git_commit"=>audit["git_commit"],"runtime_seconds"=>time()-started,"correlation_length"=>nothing,
            "sector_interpretation"=>"Exact full-spin variational calculation in an empirically low-energy loop sector; compare to unrestricted runs, not a proof that this sector is globally optimal at all lengths.")
        atomic_json(stem*".json",completed_finite_checkpoint(psi,result,cap,"projected0"))
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    for L in parse.(Int,ARGS);run_projected_zigzag(L,512);end
end
