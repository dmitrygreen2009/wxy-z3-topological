include("cylinders.jl")
include("continue.jl")
include("charge_sectors.jl")
include("observables.jl")
include("cylinder_loops.jl")
# Execute after resume_finite.jl, using its saved states rather than restarting.
# Zigzag width one has parallel-edge identifications; add a fourth width so fits can
# also omit that narrowest quotient.
for family in ["zigzag","armchair"]
    run_cylinder(family,4,4,256;ordering="star")
end
for family in ["zigzag","armchair"]
    continue_cylinder("results/$(family)_L4_w3_chi256_star.jls",512)
    run_cylinder(family,6,3,256;ordering="star")
end
for family in ["zigzag","armchair"]
    for w=1:2
        source="results/$(family)_L4_w$(w)_chi256_star.jls"
        continue_cylinder(source,512)
    end
    run_cylinder(family,6,1,256;ordering="star")
    run_cylinder(family,10,1,256;ordering="star")
    continue_cylinder("results/$(family)_L10_w1_chi256_star.jls",512)
    charge_sectors("results/$(family)_L10_w1_chi512_star.jls";chi=512)
    run_cylinder(family,8,2,256;ordering="star")
    for path in ["results/$(family)_L4_w2_chi512_star.jls","results/$(family)_L10_w1_chi512_star.jls","results/$(family)_L8_w2_chi256_star.jls"]
        measure_file(path)
    end
    loop_expectation("results/$(family)_L4_w2_chi512_star.jls")
end
