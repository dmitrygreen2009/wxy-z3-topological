include("infinite.jl")
function continue_infinite(path;maxchi=nothing,passes=2,maxiter=40)
meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
psi=load_state(path);family=meta["family"];w=meta["width"]
ordering=get(meta,"infinite_ordering","matter_first");tag=get(meta,"measurement_tag","")
ss=siteinds(only,psi);H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering),ss)
Random.seed!(7103)
target=maxchi===nothing ? max(meta["cap"],128) : maxchi
@assert target>=meta["cap"] && passes>=1 && maxiter>=1
caps=sort(unique(vcat(meta["cap"],[c for c in [32,64,128,256,512] if meta["cap"]<c<target],target)))
for cap in caps
    for pass=1:passes
        estimate_memory(psi,cap;label="$(family) infinite w$(w) expansion")
        psi=subspace_expansion(psi,H;cutoff=1e-10,maxdim=cap)
        psi=audited_vumps(H,psi;family,w,cap,maxiter,tol=1e-7,solver_tol=x->1e-10,ordering,tag)
        measure_infinite(psi,H,family,w,cap,100+10cap+pass;tag,ordering)
    end
end

end
if abspath(PROGRAM_FILE)==@__FILE__;continue_infinite(ARGS[1]);end
