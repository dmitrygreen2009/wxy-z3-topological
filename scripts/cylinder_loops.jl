isdefined(Main,:cylinder) || include("../src/model.jl")
using Serialization

function loop_expectation(path)
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    family=meta["family"];L=meta["length"];w=meta["width"]
    @assert w>=2
    lat=cylinder(family,L,w;ordering=get(meta,"ordering","axial"))
    psi=load_state(path);sites=siteinds(psi);ip=invperm(lat.order)
    wrap(x,y)=family=="zigzag" ? (x,mod(y,w)) : (x-y,mod(y,w))
    coord(t,j)=family=="zigzag" ? (t,j) : (t+j,j)
    lookup=Dict(v=>k for (k,v) in enumerate(lat.verts))
    results=[]
    for j=0:w-1
        x,y=coord(L÷2,j)
        descriptors=[(x,y,1),(x,y+1,3),(x,y+1,2),(x-1,y+1,1),(x-1,y+1,3),(x,y,2)]
        p=zeros(Int,lat.n-3lat.nv)
        for (k,(ax,ay,leg)) in enumerate(descriptors)
            t,jj=wrap(ax,ay);v=lookup[(t,jj,0)]
            e=lat.legs[v][leg];p[e]=mod(p[e]+(isodd(k) ? 1 : 2),3)
        end
        @assert all(mod(sum(p[e] for e in es),3)==0 for es in lat.legs)
        gates=ITensor[]
        for v=1:lat.nv
            c=p[lat.legs[v][1]];q=mod(p[lat.legs[v][2]]-c,3)
            (c==0 && q==0) && continue
            u=zeros(ComplexF64,8,8)
            for old=0:7
                new=sum(((old>>a)&1)<<mod(a+q,3) for a=0:2)
                u[new+1,old+1]=cis(2pi*c*(3-count_ones(old))/3)
            end
            ss=[sites[ip[3(v-1)+a]] for a=1:3]
            push!(gates,ITensor(reshape(u,2,2,2,2,2,2),prime.(ss)...,dag.(ss)...))
        end
        for e in eachindex(p)
            p[e]==0 && continue
            s=sites[ip[3lat.nv+e]]
            push!(gates,ITensor(ComplexF64[cis(2pi*p[e]/3) 0;0 1],s',dag(s)))
        end
        rotated=apply(gates,psi;cutoff=1e-13,maxdim=4maxlinkdim(psi))
        z=inner(psi,rotated)/inner(psi,psi)
        push!(results,Dict("j"=>j,"edge_exponents"=>p,"expectation_real"=>real(z),"expectation_imag"=>imag(z),
            "magnitude"=>abs(z),"rotated_norm"=>norm(rotated),"maxlinkdim_after_apply"=>maxlinkdim(rotated)))
    end
    open(replace(path,".jls"=>"_loops.json"),"w") do io;JSON3.write(io,Dict("source"=>path,"plaquette_symmetry_expectations"=>results));end
end
if abspath(PROGRAM_FILE)==@__FILE__;loop_expectation(ARGS[1]);end
