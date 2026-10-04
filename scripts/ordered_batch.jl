include("cylinders.jl")
for family=["zigzag","armchair"]
    run_cylinder(family,6,2,128;ordering="star")
end
