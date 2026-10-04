include("infinite.jl")
for family=["zigzag","armchair"]
    run_infinite(family,1,64)
end
