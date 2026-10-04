include("infinite.jl")
for path in ARGS
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    psi=load_state(path);family=meta["family"];w=meta["width"]
    ss=siteinds(only,psi)
    H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering=get(meta,"infinite_ordering","matter_first")),ss)
    measure_infinite(psi,H,family,w,meta["cap"],meta["iteration"];tag=get(meta,"measurement_tag",""),ordering=get(meta,"infinite_ordering","matter_first"))
end
