isdefined(Main,:cylinder) || include("../src/model.jl")
using Serialization
function charge_sectors(path;chi=256,deltas=[-1,1,-2,2,-3,3],tag="")
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    ground=load_state(path);family=meta["family"];L=meta["length"];w=meta["width"]
    ordering=get(meta,"ordering","axial");seed=get(meta,"seed",7103)
    lat=cylinder(family,L,w;ordering)
    sites=siteinds(ground);H=MPO(lat.os,sites)
    E0=real(inner(ground',H,ground));records=[]
    reference_nup=get(meta,"nup",cld(lat.n,2))
    for delta in deltas
        0<=reference_nup+delta<=lat.n || continue
        started=time();operation=delta>0 ? "S+" : "S-"
        sector_seed=seed+100+delta;Random.seed!(sector_seed)
        target_nup=reference_nup+delta
        partner=findlast(r->r["nup"]==lat.n-target_nup,records)
        if target_nup==lat.n-reference_nup
            psi=particle_hole_mps(ground)
            initialization="Exact particle-hole partner of reference $path"
        elseif partner!==nothing
            parent=records[partner]["checkpoint_file"]
            psi=particle_hole_mps(load_state(parent))
            initialization="Exact particle-hole partner of previously optimized $parent"
        else
            psi=deepcopy(ground);v=fld(lat.nv,2)
            targets=sort(collect(1:3lat.nv);by=m->abs(fld(m-1,3)+1-v))
            for a=1:abs(delta)
                k=findfirst(==(targets[a]),lat.order)
                orthogonalize!(psi,k)
                psi[k]=noprime(op(operation,sites[k])*psi[k])
            end
            initialization="Apply $operation at $(abs(delta)) matter sites to $path"
            if norm(psi)<=1e-12
                pattern=vcat(fill("Up",target_nup),fill("Dn",lat.n-target_nup))
                shuffle!(pattern);psi=random_mps(ComplexF64,sites,pattern;linkdims=8)
                initialization="Random QN MPS fallback after zero-norm matter-operator initialization"
            end
        end
        normalize!(psi)
        @assert val(flux(psi),"Sz")==2(reference_nup+delta)-lat.n
        estimate_memory(psi,chi;label="$(family) L$(L) w$(w) delta Nup=$(delta)")
        audit=run_provenance(;seed=sector_seed,solver="ITensorMPS charge-sector DMRG",settings=Dict("sweeps"=>12,"chi"=>chi,"cutoff"=>1e-11,"krylovdim"=>12,"delta_nup"=>delta),
            initialization,conserved_quantum_numbers=["total N_up=$(reference_nup+delta)"])
        obs=finite_checkpoint_observer(lat;family,L,w,cap=chi,seed=sector_seed,stage="charge_delta$(delta)",ordering,audit,cutoff=1e-11,noise=[1e-7,1e-8,0],nup=reference_nup+delta)
        e,psi=dmrg(H,psi;nsweeps=12,maxdim=chi,cutoff=1e-11,noise=[1e-7,1e-8,0],
            eigsolve_krylovdim=12,observer=obs,outputlevel=1)
        S,p=entropy_at(psi,lat.cut)
        record=Dict("family"=>family,"length"=>L,"width"=>w,"ordering"=>ordering,"physical_circumference"=>lat.circumference,"cap"=>chi,"solver_settings"=>audit["solver_settings"],"seed"=>sector_seed,"delta_nup"=>delta,"nup"=>reference_nup+delta,"energy"=>e,"energy_relative_to_reference"=>e-E0,
            "entropy"=>S,"sweep_energies"=>energies(obs),"truncation_errors"=>truncerrors(obs),"audit"=>audit,"runtime_seconds"=>time()-started,
            "physical_spins"=>lat.n,"bond_dimension"=>maxlinkdim(psi),"schmidt_probabilities"=>p)
        checkpoint="results/checkpoints/$(family)_w$(w)_L$(L)_N$(lat.n)_chi$(chi)_Nup$(record["nup"])_seed$(sector_seed)_v$(RUN_FORMAT_VERSION)_$(ordering)_charge_complete.jls"
        push!(records,save_checkpoint(checkpoint,psi,merge(copy(record),Dict("kind"=>"finite"))))
        result=Dict("reference_file"=>path,"reference_energy"=>E0,"chi"=>chi,"geometry"=>Dict("family"=>family,"length"=>L,"width"=>w,"physical_spins"=>lat.n,"physical_circumference"=>lat.circumference,"reference_nup"=>reference_nup),"records"=>records,
            "interpretation"=>"Finite variational sector comparison; Unexamined number and CGS sectors remain possible lower-energy competitors. Charge-three matter operator is CGS invariant.")
        atomic_json(replace(path,".jls"=>"_charge_sectors$(tag).json"),result)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;charge_sectors(ARGS[1]);end
