isdefined(Main,:cylinder) || include("../src/model.jl")
using Serialization
function measure_file(path)
    psi=load_state(path)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String))
    family=String(meta.family);L=Int(meta.length);w=Int(meta.width)
    lat=cylinder(family,L,w;ordering=haskey(meta,:ordering) ? String(meta.ordering) : "axial");invorder=invperm(lat.order)
    # Same matter species and sublattice, j=0, across the axial direction.
    vertices=[findfirst(==((t,0,0)),lat.verts) for t=0:L-1]
    selected=[invorder[3(v-1)+1] for v in vertices]
    sz=expect(psi,"Sz";sites=selected)
    zz=correlation_matrix(psi,"Sz","Sz";sites=selected)
    pm=correlation_matrix(psi,"S+","S-";sites=selected)
    connected=real.(zz).-real.(sz)*transpose(real.(sz))
    all_matter=[invorder[3(v-1)+a] for v in vertices for a=1:3]
    matter_mean=real.(expect(psi,"Sz";sites=all_matter))
    matter_zz=real.(correlation_matrix(psi,"Sz","Sz";sites=all_matter))
    star_mean=[sum(matter_mean[3(t-1)+1:3t]) for t=1:L]
    star_connected=[sum(matter_zz[3(t-1)+1:3t,3(u-1)+1:3u])-star_mean[t]*star_mean[u] for t=1:L,u=1:L]
    gauge_sites=[invorder[3lat.nv+lat.legs[v][1]] for v in vertices]
    gauge_mean=real.(expect(psi,"Sz";sites=gauge_sites))
    gauge_zz=real.(correlation_matrix(psi,"Sz","Sz";sites=gauge_sites))
    gauge_connected=gauge_zz-gauge_mean*transpose(gauge_mean)
    entropies=[entropy_at(psi,b)[1] for b=1:length(psi)-1]
    r=Dict("source"=>path,"matter_sample_sites"=>selected,"sz"=>real.(sz),
        "connected_sz"=>[collect(row) for row in eachrow(connected)],
        "splus_sminus_real"=>[collect(row) for row in eachrow(real.(pm))],
        "splus_sminus_imag"=>[collect(row) for row in eachrow(imag.(pm))],
        "entropy_all_bonds"=>entropies,
        "star_total_matter_sz"=>star_mean,
        "connected_star_total_matter_sz"=>[collect(row) for row in eachrow(star_connected)],
        "gauge_sample_sites"=>gauge_sites,"gauge_sz"=>gauge_mean,
        "connected_gauge_sz"=>[collect(row) for row in eachrow(gauge_connected)],
        "selection_rule_note"=>"Single-species matter correlators need not be invariant under CGS; star total matter Sz and gauge Sz are invariant.")
    open(replace(path,".jls"=>"_observables.json"),"w") do io;JSON3.write(io,r);end
end
if abspath(PROGRAM_FILE)==@__FILE__;for path in ARGS;measure_file(path);end;end
