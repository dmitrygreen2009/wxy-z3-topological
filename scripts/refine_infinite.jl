isdefined(Main,:measure_infinite) || include("infinite.jl")
function refine_infinite(path;rounds=3)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    family=meta["family"];w=meta["width"];cap=meta["cap"]
    psi=load_state(path);ss=siteinds(only,psi)
    H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering=get(meta,"infinite_ordering","matter_first")),ss)
    Random.seed!(7103)
    tag=get(meta,"measurement_tag","")*"_refined"
    ordering=get(meta,"infinite_ordering","matter_first")
    for k=1:rounds
        psi=audited_vumps(H,psi;family,w,cap,tol=1e-7,maxiter=20,solver_tol=x->1e-10,ordering,tag)
        measure_infinite(psi,H,family,w,cap,meta["iteration"]+k;tag,ordering)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;refine_infinite(ARGS[1]);end
