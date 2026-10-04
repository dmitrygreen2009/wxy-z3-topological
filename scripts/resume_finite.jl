include("cylinders.jl")
include("continue.jl")
for family=["zigzag","armchair"]
    continue_cylinder("results/$(family)_L6_w2_chi128_star.jls",256)
end
for w=1:3,family=["zigzag","armchair"]
    run_cylinder(family,4,w,256;ordering="star")
end
