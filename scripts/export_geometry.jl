include("../src/model.jl")
include("../src/infinite_model.jl")
mkpath("geometry")
for (name,legs) in [("single_star",[[1,2,3]]),("shared_edge_pair",[[1,2,3],[1,4,5]])]
    open("geometry/$(name).json","w") do io;JSON3.write(io,explicit_spin_graph(name,legs;vertices=[Dict("id"=>v,"x"=>0,"y"=>0,"sublattice"=>v==1 ? "A" : "B") for v=1:length(legs)],conventions=Dict("periodic_identifications"=>nothing,"B_endpoint_uses_same_leg"=>true,"A_leg_displacements_to_B"=>[[0,0],[-1,0],[0,-1]])));end
end
legs=[[1,2,3],[1,5,3],[4,5,6],[4,2,6]]
vertices=[Dict("id"=>v,"x"=>fld(v-1,2),"y"=>0,"sublattice"=>isodd(v) ? "A" : "B") for v=1:4]
convention=Dict("periodic_translations_xy"=>[[2,0],[0,1]],"A_leg_displacements_to_B"=>[[0,0],[-1,0],[0,-1]],
    "B_endpoint_uses_same_leg"=>true,"parallel_edges_distinct"=>true)
open("geometry/torus_four_star.json","w") do io;JSON3.write(io,explicit_spin_graph("torus_four_star",legs;vertices,conventions=convention));end
for family=["zigzag","armchair"],L=[2,4,6,8,10],w=1:4
    export_cylinder_geometry(family,L,w)
end
for family=["zigzag","armchair"],w=1:4,ordering=["matter_first","star"]
    data=Dict("family"=>family,"width"=>w,"ordering"=>ordering,"cell_slices"=>2,"cell_spins"=>18w,
        "translation_site_offset"=>18w,"pairs"=>infinite_pairs(family,w;ordering),
        "parent_cylinder_geometry"=>"$(family)_L6_w$(w)_star.json",
        "periodic_site_rule"=>"site k+n is translated one two-slice cell; do not identify different physical sites within the cell")
    open("geometry/infinite_$(family)_w$(w)_$(ordering).json","w") do io;JSON3.write(io,data);end
end
println("Explicit benchmark and cylinder geometries exported")
