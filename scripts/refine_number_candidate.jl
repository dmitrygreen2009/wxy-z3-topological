include("continue.jl")
# Preserve the original variational scan and expose its lowest candidate through
# a normal checksum-aware finite-result alias, retaining its actual N_up.
neutral="results/armchair_L4_w1_chi256_star.jls"
meta=JSON3.read(read(replace(neutral,".jls"=>".json"),String),Dict{String,Any})
scan=JSON3.read(read("results/armchair_L4_w1_chi256_star_charge_sectors_global_audit.json",String),Dict{String,Any})
candidate=scan["records"][argmin([r["energy"] for r in scan["records"]])]
psi=load_state(candidate["checkpoint_file"])
@assert val(flux(psi),"Sz")==2candidate["nup"]-candidate["physical_spins"]
meta["nup"]=candidate["nup"];meta["explicit_number_sector"]=true;meta["seed"]=candidate["seed"]
meta["checkpoint_file"]=candidate["checkpoint_file"]
meta["records"]=[Dict("cap"=>candidate["cap"],"energy"=>candidate["energy"],"entropy"=>candidate["entropy"],
    "maxlinkdim"=>candidate["bond_dimension"],"schmidt_probabilities"=>candidate["schmidt_probabilities"],
    "source_scan"=>"results/armchair_L4_w1_chi256_star_charge_sectors_global_audit.json")]
meta["sector_candidate_audit"]=candidate["audit"]
alias="results/armchair_L4_w1_chi256_star_number_candidate_Nup$(candidate["nup"]).jls"
atomic_json(replace(alias,".jls"=>".json"),meta)
continue_cylinder(neutral,512;cutoff=1e-12,noise_first=[1e-8,1e-9,0])
continue_cylinder(alias,512;cutoff=1e-12,noise_first=[1e-8,1e-9,0])
