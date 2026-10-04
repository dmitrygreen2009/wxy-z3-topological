include("infinite.jl")
include("observables.jl")
include("cylinder_loops.jl")
path="results/infinite_zigzag_w1_chi32.jls"
meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
psi=load_state(path);ss=siteinds(only,psi)
H=InfiniteSum{MPO}(infinite_opsum("zigzag",1),ss)
measure_infinite(psi,H,"zigzag",1,32,meta["iteration"])
for path in ["results/armchair_L2_w1_chi256_star.jls","results/zigzag_L6_w1_chi256.jls","results/zigzag_L6_w2_chi256_star.jls"]
    measure_file(path)
end
loop_expectation("results/zigzag_L6_w2_chi256_star.jls")
