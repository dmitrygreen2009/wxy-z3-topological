include("cylinders.jl")
for family=["zigzag","armchair"],seed=[7104,7105]
    run_cylinder(family,4,1,256;ordering="star",seed)
end
