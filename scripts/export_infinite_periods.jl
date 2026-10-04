include("../src/model.jl")
include("../src/infinite_model.jl")
for family in ["zigzag","armchair"],w in [1,2,3],slices in [3,6],ordering in ["star","matter_first"]
 n=9w*slices
 r=Dict("family"=>family,"width"=>w,"ordering"=>ordering,"cell_slices"=>slices,"cell_spins"=>n,
 "translation_site_offset"=>n,"pairs"=>infinite_pairs(family,w;ordering,cell_slices=slices),
 "primitive_geometry"=>"primitive_honeycomb_cell.json","periodic_site_rule"=>"site k+n is translated by cell_slices spatial slices; endpoint legs are unchanged",
 "spatial_slice_cut_bonds"=>collect(9w:9w:n),"primitive_vertices_per_cell"=>2w*slices,
 "physical_circumference"=>family=="zigzag" ? sqrt(3)*w : 3.0*w,
 "physical_axial_cell_length"=>slices*(family=="zigzag" ? 1.5 : sqrt(3)/2),
 "A_B_W_convention"=>"Same W at both endpoints, columns labeled by physical leg i")
 open("geometry/infinite_$(family)_w$(w)_cell$(slices)_$(ordering).json","w") do io;JSON3.write(io,r);end
end
