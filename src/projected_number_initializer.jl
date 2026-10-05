# Exact number-only randomization followed by the validated P_q B† bridge.
# This avoids the extra occupancy constraints of adjacent two-site gates with
# unequal winding weights. All tensor operations use ITensorMPS.
isdefined(@__MODULE__,:bridge_physical_checkpoint) || include("winding_checkpoint_bridge.jl")
function projected_number_initializer(lat,p,nup,charge;seed=7254,linkdims=8,maxdim=128,cutoff=1e-13)
    Random.seed!(seed)
    physical_sites=siteinds("S=1/2",lat.n;conserve_qns=true)
    state=shuffle(vcat(fill("Up",nup),fill("Dn",lat.n-nup)))
    physical=random_mps(Float64,physical_sites,state;linkdims)
    converted,record=bridge_physical_checkpoint(physical,lat,p;charge,maxdim,cutoff)
    @assert record["projection_weight"]>1e-12
    @assert val(flux(converted),"Sz")==2nup-lat.n && mod(val(flux(converted),"Winding"),3)==charge
    record["strategy"]="Physical number-only random MPS followed by exact microscopic P_q B†"
    record["physical_initial_product_state"]=state
    record["physical_random_linkdims"]=linkdims
    record["physical_random_scalar_type"]="Float64"
    record["converted_scalar_type"]=string(eltype(converted[1]))
    record["seed"]=seed
    converted,record
end
