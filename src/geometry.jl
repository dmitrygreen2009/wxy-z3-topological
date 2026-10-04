function explicit_spin_graph(name,legs;vertices=nothing,conventions=Dict(),order=nothing,cut=nothing)
    nv=length(legs);ng=maximum(vcat(legs...));n=3nv+ng
    vertices===nothing && (vertices=[Dict("id"=>v) for v=1:nv])
    matter=[Dict("physical_spin_id"=>3(v-1)+a,"vertex_id"=>v,"species_a"=>a,"spin"=>"1/2") for v=1:nv for a=1:3]
    edges=[Dict("edge_id"=>e,"physical_spin_id"=>3nv+e,"spin"=>"1/2",
        "endpoints"=>[Dict("vertex_id"=>v,"local_leg_i"=>i,"W_convention"=>"W[a,i]") for v=1:nv for i=1:3 if legs[v][i]==e],
        "dangling"=>count(==(e),vcat(legs...))==1) for e=1:ng]
    for edge in edges
        endpoints=edge["endpoints"]
        av=findfirst(ep->get(vertices[ep["vertex_id"]],"sublattice","")=="A",endpoints)
        bv=findfirst(ep->get(vertices[ep["vertex_id"]],"sublattice","")=="B",endpoints)
        (av===nothing && bv===nothing) && continue
        ep=av===nothing ? endpoints[bv] : endpoints[av];v=vertices[ep["vertex_id"]]
        x=get(v,"bravais_x",get(v,"x",0));y=get(v,"bravais_y",get(v,"y",0));leg=ep["local_leg_i"]
        ax,ay=av===nothing ? (x+(leg==2 ? 1 : 0),y+(leg==3 ? 1 : 0)) : (x,y)
        bx,by=ax-(leg==2 ? 1 : 0),ay-(leg==3 ? 1 : 0)
        edge["owner_A_unwrapped_bravais"]=[ax,ay]
        edge["target_B_unwrapped_bravais"]=[bx,by]
        edge["physical_midpoint_xy"]=[sqrt(3)*(ax-ay+bx-by)/4,0.75*(ax+ay+bx+by)+0.5]
        if bv!==nothing
            vB=vertices[endpoints[bv]["vertex_id"]]
            edge["periodic_translation_to_B_representative"]=[get(vB,"bravais_x",get(vB,"x",0))-bx,get(vB,"bravais_y",get(vB,"y",0))-by]
        end
    end
    Dict("name"=>name,"index_base"=>1,"physical_spins"=>n,"vertices"=>vertices,"matter_spins"=>matter,
        "shared_gauge_edges"=>edges,"endpoint_leg_table"=>legs,"conventions"=>conventions,
        "mps_order_physical_spin_ids"=>order,"spatial_cut_mps_bond"=>cut)
end
function export_cylinder_geometry(family,L,w;ordering="star")
    lat=cylinder(family,L,w;ordering)
    vertices=[Dict("id"=>v,"t"=>t,"j"=>j,"sublattice"=>s==0 ? "A" : "B",
        "bravais_x"=>family=="zigzag" ? t : t+j,"bravais_y"=>j) for (v,(t,j,s)) in enumerate(lat.verts)]
    convention=Dict("family"=>family,"length_slices"=>L,"width_cells"=>w,"circumference"=>lat.circumference,
        "bravais_a1"=>[sqrt(3)/2,1.5],"bravais_a2"=>[-sqrt(3)/2,1.5],"B_minus_A"=>[0,1],
        "periodic_translation_xy"=>family=="zigzag" ? [0,w] : [w,w],"open_axis"=>true,
        "A_leg_displacements_to_B"=>[[0,0],[-1,0],[0,-1]],"B_endpoint_uses_same_leg"=>true,
        "dangling_boundary_gauge_spins_retained"=>true,"ordering"=>ordering,
        "spatial_partition_rule"=>"Half-open physical axial coordinate u < L/2, in spatial-slice units",
        "physical_axial_matter_offsets_A_B"=>family=="zigzag" ? [0,1/3] : [0,0],
        "physical_axial_gauge_offsets_from_owner_A_by_leg"=>family=="zigzag" ? [1/6,-1/3,1/6] : [0,-1/2,1/2],
        "ordering_partition_proxy"=>"Matter key t+0.01a; gauge key (owner_A_t+target_B_t)/2+0.1. For integer L these keys give the same half-open partition as the physical coordinates; they are not physical length coordinates.")
    graph=explicit_spin_graph("$(family)_L$(L)_w$(w)",lat.legs;vertices,conventions=convention,order=lat.order,cut=lat.cut)
    open("geometry/$(family)_L$(L)_w$(w)_$(ordering).json","w") do io;JSON3.write(io,graph);end
end
