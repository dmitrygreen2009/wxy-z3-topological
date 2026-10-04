include("cylinders.jl")
for family=["zigzag","armchair"]
    run_cylinder(family,10,1,512)
    run_cylinder(family,4,3,256)
end
