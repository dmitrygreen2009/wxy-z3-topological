include("cylinders.jl")
for w=1:3, family=["zigzag","armchair"]
    run_cylinder(family,4,w,256;ordering="star")
end
