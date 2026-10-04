isdefined(Main,:cylinder) || include("../src/model.jl")
using Serialization
function charge_sectors(path;chi=256,deltas=[-1,1,-3,3])
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    ground=load_state(path);family=meta["family"];L=meta["length"];w=meta["width"]
    ordering=get(meta,"ordering","axial");seed=get(meta,"seed",7103)
    lat=cylinder(family,L,w;ordering)
    sites=siteinds(ground);H=MPO(lat.os,sites)
    E0=real(inner(ground',H,ground));records=[]
    for delta in deltas
        started=time();operation=delta>0 ? "S+" : "S-"
        psi=deepcopy(ground);v=fld(lat.nv,2)
        for a=1:abs(delta)
            k=findfirst(==(3(v-1)+a),lat.order)
            orthogonalize!(psi,k)
            psi[k]=noprime(op(operation,sites[k])*psi[k])
        end
        norm(psi)>1e-12 || error("Zero-norm charge initialization for delta=$delta; choose another matter star")
        normalize!(psi)
        estimate_memory(psi,chi;label="$(family) L$(L) w$(w) delta Nup=$(delta)")
        audit=run_provenance(;seed,solver="ITensorMPS charge-sector DMRG",settings=Dict("sweeps"=>12,"chi"=>chi,"cutoff"=>1e-11,"krylovdim"=>12,"delta_nup"=>delta),
            initialization="Apply $operation at $(abs(delta)) matter sites to $path",conserved_quantum_numbers=["total N_up=$(cld(lat.n,2)+delta)"])
        obs=finite_checkpoint_observer(lat;family,L,w,cap=chi,seed,stage="charge_delta$(delta)",ordering,audit,cutoff=1e-11,noise=[1e-7,1e-8,0],nup=cld(lat.n,2)+delta)
        e,psi=dmrg(H,psi;nsweeps=12,maxdim=chi,cutoff=1e-11,noise=[1e-7,1e-8,0],
            eigsolve_krylovdim=12,observer=obs,outputlevel=1)
        S,p=entropy_at(psi,lat.cut)
        record=Dict("delta_nup"=>delta,"nup"=>cld(lat.n,2)+delta,"energy"=>e,"energy_relative_to_reference"=>e-E0,
            "entropy"=>S,"sweep_energies"=>energies(obs),"truncation_errors"=>truncerrors(obs),"audit"=>audit,"runtime_seconds"=>time()-started,
            "physical_spins"=>lat.n,"bond_dimension"=>maxlinkdim(psi),"schmidt_probabilities"=>p)
        checkpoint="results/checkpoints/$(family)_w$(w)_L$(L)_N$(lat.n)_chi$(chi)_Nup$(record["nup"])_seed$(seed)_v$(RUN_FORMAT_VERSION)_$(ordering)_charge_complete.jls"
        push!(records,save_checkpoint(checkpoint,psi,merge(copy(record),Dict("kind"=>"finite"))))
        result=Dict("reference_file"=>path,"reference_energy"=>E0,"chi"=>chi,"records"=>records,
            "interpretation"=>"Finite variational sector comparison; no unexamined number or flux sectors excluded. Charge-three matter operator is CGS invariant.")
        atomic_json(replace(path,".jls"=>"_charge_sectors.json"),result)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;charge_sectors(ARGS[1]);end
