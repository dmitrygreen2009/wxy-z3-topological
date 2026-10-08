# Measurement-only 2D trial: product of reduced rectangular cylinder patches.
include("../src/model.jl")
function tiled_trial(path,measurement)
 started=time();meta=JSON3.read(read(path,String),Dict{String,Any});saved=JSON3.read(read(measurement,String),Dict{String,Any})
 f=meta["family"];L=meta["length"];w=meta["width"];lat=cylinder(f,L,w;ordering="star")
 @assert meta["mps_order"]==lat.order && saved["physical_spins"]==lat.n
 checkpoint=saved["source_checkpoint"];digest=open(io->bytes2hex(sha256(io)),checkpoint);@assert digest==saved["source_sha256"]
 psi=load_state(checkpoint);@assert length(psi)==lat.n && abs(norm(psi)-1)<1e-6
 @assert (val(flux(psi),"Sz")+lat.n)/2==saved["nup"]
 geo=JSON3.read(read("geometry/$(f)_L$(L)_w$(w)_star.json",String),Dict{String,Any})
 verts=Dict(v["id"]=>v for v in geo["vertices"]);cells=Set((v["bravais_x"],v["bravais_y"]) for v in values(verts))
 owners=Dict{Int,Int}();targets=Dict{Int,Tuple{Int,Int}}();keep=Int[];shifts=[(0,0),(-1,0),(0,-1)]
 for e in geo["shared_gauge_edges"]
  ends=[x for x in e["endpoints"] if verts[x["vertex_id"]]["sublattice"]=="A"]
  isempty(ends) && continue
  @assert length(ends)==1
  a=only(ends);v=verts[a["vertex_id"]];i=a["local_leg_i"];edge=e["edge_id"]
  owners[edge]=v["id"];targets[edge]=(v["bravais_x"]+shifts[i][1],v["bravais_y"]+shifts[i][2]);push!(keep,e["physical_spin_id"])
 end
 @assert length(keep)==3L*w && length(unique(keep))==length(keep)
 append!(keep,1:3lat.nv);sort!(keep);@assert length(keep)==9L*w
 os=OpSum();ip=invperm(lat.order);retained=[];removed=[]
 for v=1:lat.nv,i=1:3
  e=lat.legs[v][i];inside=haskey(owners,e) && (verts[v]["sublattice"]=="A" || targets[e] in cells)
  if inside
   if verts[v]["sublattice"]=="B";@assert targets[e]==(verts[v]["bravais_x"],verts[v]["bravais_y"]);end
   push!(retained,(v,i,e))
   for a=1:3
    m=ip[3(v-1)+a];g=ip[3lat.nv+e]
    os += -W[a,i],"S-",m,"S+",g
    os += -conj(W[a,i]),"S+",m,"S-",g
   end
  else
   push!(removed,(v,i,e))
  end
 end
 energy=inner(psi',MPO(os,siteinds(psi)),psi)/inner(psi,psi);@assert abs(imag(energy))<1e-10
 @assert open(io->bytes2hex(sha256(io)),checkpoint)==digest
 density=saved["physical_site_density_profile"];N=sum(density[j] for j in keep)
 result=Dict("source_result"=>path,"source_measurement"=>measurement,"source_checkpoint"=>checkpoint,"source_sha256"=>digest,"family"=>f,"length"=>L,"width"=>w,
 "primitive_cells_per_patch"=>L*w,"patch_physical_spin_count"=>length(keep),"original_physical_spins"=>lat.n,"kept_physical_spin_ids"=>keep,
 "retained_star_leg_terms"=>retained,"removed_star_leg_terms"=>removed,"patch_energy"=>real(energy),"energy_imaginary_error"=>abs(imag(energy)),"bulk_trial_energy_per_primitive_cell"=>real(energy)/(L*w),
 "patch_nup_expectation"=>N,"patch_sz_expectation"=>N-length(keep)/2,"patch_filling_fraction"=>N/length(keep),"original_nup"=>saved["nup"],"original_sz"=>saved["total_sz"],
 "state_construction"=>"Keep all patch matter and exactly the gauges owned by patch A vertices; trace out extra outside-owner gauges. Tile the resulting reduced U1-invariant density operator as a product over the infinite plane. Gauge edges are assigned once to A owner; remove original periodic-wrap B couplings and external-owner B couplings. Cross-patch hopping expectations vanish exactly because every retained one-spin ladder expectation is zero by U1.",
 "claim"=>"Valid 2D global-ground energy upper bound from a possibly unconverged trial, not evidence for its phase or its pure ground filling. No Hamiltonian is changed: omitted interactions exist in the tiled physical Hamiltonian and have zero expectation.",
 "mixed_trial"=>true,"optimization_performed"=>false,"runtime_seconds"=>time()-started,"git_commit"=>LAUNCH_REVISION)
 atomic_json("results/physics48_tiled_trial.json",result);println("2D tiled trial per primitive cell ",real(energy)/(L*w)," rho=",N/length(keep)," kept spin count ",length(keep));flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__
 tiled_trial("results/zigzag_L6_w3_chi256_star.json","results/physics48_zigzag_L6_w3_chi256.json")
end
