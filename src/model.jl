using ITensors, ITensorMPS, LinearAlgebra, Random, JSON3
BLAS.set_num_threads(1)
if Threads.nthreads()>1
    ITensors.NDTensors.Strided.disable_threads()
    ITensors.enable_threaded_blocksparse()
end
const W = ComplexF64[1 1 1; 1 cis(2pi/3) cis(-2pi/3); 1 cis(-2pi/3) cis(2pi/3)] / sqrt(3)
struct LoggedObserver <: AbstractObserver
    data::DMRGObserver
    checkpoint_callback::Union{Nothing,Function}
end
LoggedObserver(;checkpoint_callback=nothing)=LoggedObserver(DMRGObserver(),checkpoint_callback)
function ITensorMPS.measure!(o::LoggedObserver;kwargs...)
    ITensorMPS.measure!(o.data;kwargs...)
    if get(kwargs,:sweep_is_done,false)
        flush(stdout)
        o.checkpoint_callback===nothing || o.checkpoint_callback(kwargs)
    end
end
function ITensorMPS.checkdone!(o::LoggedObserver;kwargs...)
    flush(stdout)
    ITensorMPS.checkdone!(o.data;kwargs...)
end
ITensorMPS.energies(o::LoggedObserver)=energies(o.data)
ITensorMPS.truncerrors(o::LoggedObserver)=truncerrors(o.data)

function microscopic_mpo(legs, sites)
    os=OpSum()
    for v in eachindex(legs), a=1:3, i=1:3
        m=3(v-1)+a; g=3length(legs)+legs[v][i]
        os += -W[a,i], "S-", m, "S+", g
        os += -conj(W[a,i]), "S+", m, "S-", g
    end
    MPO(os,sites)
end

function entropy_at(psi,b)
    ITensorMPS.orthogonalize!(psi,b)
    left=b==1 ? (siteind(psi,b),) : (linkind(psi,b-1),siteind(psi,b))
    _,s,_=svd(psi[b],left)
    p=abs2.(diag(Array(s,inds(s)...)))
    p ./= sum(p)
    -sum(x>0 ? x*log(x) : 0.0 for x in p), sort(p,rev=true)
end

"""All three legs retained, with dangling physical gauge spins at open ends.
zigzag: T=(0,w); armchair: T=(w,w), in triangular Bravais coordinates.
Each leg of A(x,y) joins B(x,y), B(x-1,y), B(x,y-1).
"""
function cylinder(family,L,w;ordering="axial")
    wrap(x,y)=family=="zigzag" ? (x,mod(y,w)) : (x-y,mod(y,w))
    coord(t,j)=family=="zigzag" ? (t,j) : (t+j,j)
    verts=[(t,j,s) for t=0:L-1 for j=0:w-1 for s=0:1]
    lookup=Dict(v=>i for (i,v) in enumerate(verts))
    edges=Dict{Tuple{Int,Int,Int},Int}(); legs=[zeros(Int,3) for _ in verts]
    positions=Float64[]
    for (v,(t,j,s)) in enumerate(verts)
        x,y=coord(t,j)
        for i=1:3
            ax,ay=s==0 ? (x,y) : i==1 ? (x,y) : i==2 ? (x+1,y) : (x,y+1)
            at,aj=wrap(ax,ay); key=(at,aj,i)
            if !haskey(edges,key)
                edges[key]=length(edges)+1
                bx,by=i==1 ? (ax,ay) : i==2 ? (ax-1,ay) : (ax,ay-1)
                bt,bj=wrap(bx,by)
                push!(positions,(at+bt)/2+0.1)
            end
            legs[v][i]=edges[key]
        end
    end
    nv=length(verts); n=3nv+length(edges)
    # Individual spins sorted by axial slice; no circumference supersites.
    keys=vcat([Float64(v[1])+0.01*a for v in verts for a=1:3],positions)
    if ordering=="star"
        sortkeys=[(v[1],9v[2]+(v[3]==0 ? a : 5+a)) for v in verts for a=1:3]
        gaugekeys=Vector{Tuple{Int,Int}}(undef,length(edges))
        for ((t,j,i),e) in edges
            gaugekeys[e]=(i==2 ? t-1 : t,9j+(i==1 ? 4 : i==3 ? 5 : 9))
        end
        append!(sortkeys,gaugekeys)
        order=sortperm(1:n,by=i->(sortkeys[i],i))
    else
        order=sortperm(1:n,by=i->(keys[i],i))
    end
    invorder=invperm(order)
    os=OpSum()
    for v=1:nv,a=1:3,i=1:3
        m=invorder[3(v-1)+a]; g=invorder[3nv+legs[v][i]]
        os += -W[a,i],"S-",m,"S+",g
        os += -conj(W[a,i]),"S+",m,"S-",g
    end
    cut=count(k->k < L/2,keys)
    (;os,n,nv,legs,verts,order,cut,circumference=(family=="zigzag" ? sqrt(3)*w : 3.0*w))
end

include("audit.jl")

include("geometry.jl")
