include("../src/model.jl")
for path in ARGS
    started=time();meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    psi=load_state(path);lat=cylinder(meta["family"],meta["length"],meta["width"];ordering=get(meta,"ordering","axial"))
    H=MPO(lat.os,siteinds(psi));E=real(inner(psi',H,psi))
    variance=real(inner(H,psi,H,psi))-E^2
    @assert variance>=-1e-9 "Negative variance beyond rounding: $variance"
    S,p=entropy_at(psi,lat.cut)
    result=Dict("source"=>path,"measurement_git_commit"=>LAUNCH_REVISION,"energy"=>E,"energy_variance"=>variance,
        "eigenvector_residual_norm"=>sqrt(max(0,variance)),"full_state_entropy"=>S,"schmidt_probabilities"=>p,
        "bond_dimension"=>maxlinkdim(psi),"state_norm"=>norm(psi),"runtime_seconds"=>time()-started,
        "audit"=>run_provenance(solver="ITensorMPS exact expectation contractions",settings=Dict(),initialization=path,conserved_quantum_numbers=["total Sz"]),
        "interpretation"=>"Variance tests eigenstate quality, not global ground-state optimality; compare independent initializations and sectors.")
    atomic_json(replace(path,".jls"=>"_variance.json"),result)
    println(path," E=",E," variance=",variance," entropy=",S);flush(stdout)
end
