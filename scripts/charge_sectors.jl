isdefined(Main,:cylinder) || include("../src/model.jl")
using Serialization
function charge_sectors(path;chi=256)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    ground=load_state(path);lat=cylinder(meta["family"],meta["length"],meta["width"];ordering=get(meta,"ordering","axial"))
    sites=siteinds(ground);H=MPO(lat.os,sites)
    E0=real(inner(ground',H,ground));records=[]
    for (operation,delta) in [("S-",-1),("S+",1)]
        psi=deepcopy(ground)
        # Select a central matter site; all matter is retained in this state.
        physical=3*(fld(lat.nv,2)-1)+1
        k=findfirst(==(physical),lat.order)
        orthogonalize!(psi,k)
        psi[k]=noprime(op(operation,sites[k])*psi[k]);normalize!(psi)
        obs=LoggedObserver()
        e,psi=dmrg(H,psi;nsweeps=12,maxdim=chi,cutoff=1e-11,noise=[1e-7,1e-8,0],
            eigsolve_krylovdim=12,observer=obs,outputlevel=1)
        S,p=entropy_at(psi,lat.cut)
        push!(records,Dict("delta_nup"=>delta,"energy"=>e,"energy_relative_to_reference"=>e-E0,
            "entropy"=>S,"sweep_energies"=>energies(obs),"truncation_errors"=>truncerrors(obs)))
    end
    result=Dict("reference_file"=>path,"reference_energy"=>E0,"chi"=>chi,"records"=>records,
        "interpretation"=>"Finite variational sector comparison; excludes no unexamined number or flux sectors.")
    open(replace(path,".jls"=>"_charge_sectors.json"),"w") do io;JSON3.write(io,result);end
end
if abspath(PROGRAM_FILE)==@__FILE__;charge_sectors(ARGS[1]);end
