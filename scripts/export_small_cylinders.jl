include("../src/model.jl")
out=Dict()
for family=["zigzag","armchair"]
    lat=cylinder(family,2,1;ordering="star")
    out[family]=Dict("legs_zero_based"=>[l.-1 for l in lat.legs],"spins"=>lat.n,
        "left_physical_sites_zero_based"=>lat.order[1:lat.cut].-1,"cut"=>lat.cut)
end
open("results/small_cylinder_incidence.json","w") do io;JSON3.write(io,out);end
