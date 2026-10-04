isdefined(Main,:measure_infinite) || include("infinite.jl")
function refine_infinite(path;rounds=3)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    family=meta["family"];w=meta["width"];cap=meta["cap"]
    psi=load_state(path);ss=siteinds(only,psi)
    H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering=get(meta,"infinite_ordering","matter_first")),ss)
    for k=1:rounds
        psi=tdvp(H,psi;time_step=-Inf,tol=1e-7,maxiter=20,
            solver_tol=x->1e-10,multisite_update_alg="sequential")
        measure_infinite(psi,H,family,w,cap,meta["iteration"]+k;tag=get(meta,"measurement_tag",""),ordering=get(meta,"infinite_ordering","matter_first"))
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;refine_infinite(ARGS[1]);end
