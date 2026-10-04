isdefined(Main,:charge_sectors) || include("charge_sectors.jl")
# Compare exact number blocks using saved microscopic states; retain every point.
# Particle-hole plus matter-species 2/3 interchange gives exact mirror energies.
for (family,L,w,chi,deltas) in [("armchair",4,1,256,[-1,1,-2,2,-3,3,-4,-5,-6]),
    ("zigzag",4,1,256,collect(-1:-1:-6)),
    ("armchair",4,2,256,collect(-1:-1:-8)),
    ("zigzag",4,2,256,collect(-1:-1:-8)),
    ("armchair",6,2,256,collect(-1:-1:-8)),
    ("zigzag",6,2,256,collect(-1:-1:-8))]
    path="results/$(family)_L$(L)_w$(w)_chi256_star.jls"
    isfile(replace(path,".jls"=>".json")) || continue
    completed=Int[]
    for file in readdir("results";join=true)
        endswith(file,"_charge_sectors_global_audit.json") || continue
        data=JSON3.read(read(file,String),Dict{String,Any})
        g=get(data,"geometry",Dict())
        get(g,"family",nothing)==family && get(g,"length",nothing)==L && get(g,"width",nothing)==w && data["chi"]==chi || continue
        for record in data["records"]
            try
                state=load_state(record["checkpoint_file"])
                @assert val(flux(state),"Sz")==2record["nup"]-record["physical_spins"]
                push!(completed,record["delta_nup"])
            catch error
                println("Rejecting invalid sector checkpoint: ",sprint(showerror,error));flush(stdout)
            end
        end
    end
    pending=setdiff(deltas,completed)
    isempty(pending) && (println("Reusing complete validated sector window for $family L=$L w=$w");flush(stdout);continue)
    # Separate result suffix preserves completed points when resuming a partial window.
    tag=isempty(completed) ? "_global_audit" : "_global_audit_continuation"
    charge_sectors(path;chi,deltas=pending,tag)
end
