# Measurement only, using official ITensor contractions. Never calls DMRG/VUMPS.
include("../src/model.jl")
include("../src/cgs_operator_mpo.jl")
include("../src/finite_stationarity.jl")
function measure48(family,L,w,cap)
    started=time();path="results/$(family)_L$(L)_w$(w)_chi$(cap)_star.json"
    meta=JSON3.read(read(path,String),Dict{String,Any});payload=get(meta,"checkpoint_file",replace(path,".json"=>".jls"))
    sidecar=JSON3.read(read(replace(payload,".jls"=>".json"),String),Dict{String,Any})
    if haskey(sidecar,"checkpoint_sha256") && haskey(sidecar,"kind")
        psi,cp=load_valid_checkpoint(payload)
    else
        psi=load_state(payload);cp=sidecar
        @assert isfinite(norm(psi)) && abs(norm(psi)-1)<1e-6
    end
    digest=open(io->bytes2hex(sha256(io)),payload)
    @assert get(cp,"basis","physical_spin")=="physical_spin"
    lat=cylinder(family,L,w;ordering="star");@assert length(psi)==lat.n
    @assert meta["mps_order"]==lat.order && meta["cut"]==lat.cut
    s=siteinds(psi);ip=invperm(lat.order);sz=real.(expect(psi,"Sz"));nup=(val(flux(psi),"Sz")+lat.n)/2
    @assert abs(sum(sz)+lat.n/2-nup)<1e-8
    physical_density=(sz.+0.5)[ip];geo=JSON3.read(read("geometry/$(family)_L$(L)_w$(w)_star.json",String),Dict{String,Any})
    profiles=[]
    for t=0:L-1
        matter=[3(v-1)+a for (v,coord) in enumerate(lat.verts) if coord[1]==t for a=1:3]
        # Interior edges counted once; assign a gauge spin to the mean endpoint t.
        gauges=Int[]
        for e in geo["shared_gauge_edges"]
            ts=[lat.verts[x["vertex_id"]][1] for x in e["endpoints"]]
            sum(ts)/length(ts)==t && push!(gauges,e["physical_spin_id"])
        end
        push!(profiles,Dict("t"=>t,"matter_density_mean"=>sum(physical_density[matter])/length(matter),
         "matter_density_min"=>minimum(physical_density[matter]),"matter_density_max"=>maximum(physical_density[matter]),
         "gauge_density_mean_at_integer_endpoint_midpoint"=>isempty(gauges) ? nothing : sum(physical_density[gauges])/length(gauges)))
    end
    defs=JSON3.read(read("geometry/cgs_cycles/$(family)_L$(L)_w$(w)_star.json",String),Dict{String,Any})
    cycles=[c for c in defs["cycles"] if c["winding_number"]==1]
    # Prefer the central homologous loop, retaining its exact recorded support.
    preferred="circumference_cycle_t$(max(0,L÷2-1))"
    ix=findfirst(c->c["label"]==preferred,cycles);cycle=cycles[ix===nothing ? 1 : ix]
    p=Int.(cycle["edge_exponents"]);@assert all(mod(sum(p[e] for e in row),3)==0 for row in lat.legs)
    triplets=Tuple{Int,Int,Int}[];gauges=Tuple{Int,Int}[]
    for v=1:lat.nv
        c=mod(p[lat.legs[v][1]],3);q=mod(p[lat.legs[v][2]]-c,3)
        if c!=0 || q!=0
            positions=[ip[3(v-1)+a] for a=1:3];@assert positions==collect(first(positions):first(positions)+2)
            push!(triplets,(first(positions),c,q))
        end
    end
    for e in eachindex(p);p[e]!=0 && push!(gauges,(ip[3lat.nv+e],p[e]));end
    spec=(;triplets,gauges,first_site=1,last_site=lat.n)
    U=cgs_product_mpo(s,spec);z=inner(psi',U,psi)/inner(psi,psi)
    second=(;triplets=[(j,mod(-c,3),mod(-q,3)) for (j,c,q) in triplets],gauges=[(j,mod(-p,3)) for (j,p) in gauges],first_site=1,last_site=lat.n)
    z2=inner(psi',cgs_product_mpo(s,second),psi)/inner(psi,psi)
    @assert abs(z2-conj(z))<1e-10
    probs=[real((1+conj(cis(2pi*q/3))*z+conj(cis(2pi*q/3))^2*z2)/3) for q=0:2]
    entropy,prob=entropy_at(psi,lat.cut)
    @assert abs(entropy-meta["records"][end]["entropy"])<1e-8
    # One central projected residual, never an optimization or H^2 contraction.
    H=MPO(lat.os,s);stationarity=finite_bond_stationarity(H,psi,lat.cut)
    output=Dict("source_result"=>path,"source_checkpoint"=>payload,"source_sha256"=>digest,"family"=>family,"length"=>L,"width"=>w,"cap"=>cap,"physical_spins"=>lat.n,
      "nup"=>nup,"total_sz"=>nup-lat.n/2,"filling_fraction"=>nup/lat.n,"filling_selection"=>"Fixed by saved U1 MPS sector; not variational selection of global ground filling",
      "physical_site_density_profile"=>physical_density,"row_density_profiles"=>profiles,"spatial_cut"=>lat.cut,"full_state_entropy"=>entropy,
      "winding_cycle_label"=>cycle["label"],"winding_edge_exponents"=>p,"matter_actions"=>triplets,"gauge_actions"=>gauges,
      "winding_mean"=>[real(z),imag(z)],"winding_unitary_variance"=>1-abs2(z),"winding_charge_probabilities"=>probs,
      "winding_purity_error"=>minimum(abs(z-cis(2pi*q/3)) for q=0:2),"bond_stationarity"=>stationarity,"runtime_seconds"=>time()-started,
      "measurement_git_commit"=>LAUNCH_REVISION,"transfer_correlation_length"=>nothing,"transfer_status"=>"Finite open MPS; no arbitrary periodic wrap transfer length",
      "operator_validation"=>"Same exact local matter permutation/phase as the audited 64-state star matrices, physical shared edges once, all endpoint exponent sums zero mod3; contractions untruncated",
      "optimization_performed"=>false)
    @assert open(io->bytes2hex(sha256(io)),payload)==digest
    atomic_json("results/physics48_$(family)_L$(L)_w$(w)_chi$(cap).json",output)
    println("MEASURED ",family," L",L," w",w," cap",cap," S=",entropy," loop=",z," variance=",1-abs2(z));flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__
    cases=vcat([(family,4,w,cap) for family in ["zigzag","armchair"] for w in [2,3] for cap in [256,512]],[(family,6,w,256) for family in ["zigzag","armchair"] for w in [2,3]])
    for (family,L,w,cap) in cases
        try measure48(family,L,w,cap) catch error
            atomic_json("results/physics48_$(family)_L$(L)_w$(w)_chi$(cap)_error.json",Dict("error"=>sprint(showerror,error,catch_backtrace()),"optimization_performed"=>false));println("MEASUREMENT ERROR ",family," w",w,": ",sprint(showerror,error));flush(stdout)
        end
        GC.gc()
    end
end
