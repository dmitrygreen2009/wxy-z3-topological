include("cylinders.jl")
for family=["zigzag","armchair"]
    run_cylinder(family,6,1,256)
    run_cylinder(family,6,2,256)
end
