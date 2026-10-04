# Exact microscopic symmetry gates: physical spin rotations and matter permutations.
function cgs_cycle_gates(lat,sites,p)
    @assert length(p)==lat.n-3lat.nv
    @assert all(mod(sum(p[e] for e in es),3)==0 for es in lat.legs)
    ip=invperm(lat.order);gates=ITensor[]
    for v=1:lat.nv
        c=mod(p[lat.legs[v][1]],3);q=mod(p[lat.legs[v][2]]-c,3)
        c==0 && q==0 && continue
        u=zeros(ComplexF64,8,8)
        for old=0:7
            new=sum(((old>>a)&1)<<mod(a+q,3) for a=0:2)
            u[new+1,old+1]=cis(2pi*c*(3-count_ones(old))/3)
        end
        ss=[sites[ip[3(v-1)+a]] for a=1:3]
        push!(gates,ITensor(reshape(u,2,2,2,2,2,2),prime.(ss)...,dag.(ss)...))
    end
    for e in eachindex(p)
        mod(p[e],3)==0 && continue
        s=sites[ip[3lat.nv+e]]
        push!(gates,ITensor(ComplexF64[cis(2pi*p[e]/3) 0;0 1],s',dag(s)))
    end
    gates
end
function cgs_project(psi,lat,p;charge=0,maxdim=256,cutoff=1e-12)
    gates=cgs_cycle_gates(lat,siteinds(psi),p)
    u=apply(gates,psi;cutoff,maxdim)
    u2=apply(gates,u;cutoff,maxdim)
    z=cis(2pi*charge/3)
    projected=add(psi,conj(z)*u,conj(z)^2*u2;cutoff,maxdim)
    weight=norm(projected)^2/(9norm(psi)^2)
    weight>1e-14 || error("Requested CGS charge $charge has vanishing weight $weight")
    normalize!(projected)
    value=inner(projected,apply(gates,projected;cutoff,maxdim))/inner(projected,projected)
    projected,Dict("charge"=>charge,"projection_weight"=>weight,"symmetry_expectation_real"=>real(value),
        "symmetry_expectation_imag"=>imag(value),"purity_error"=>abs(value-z),"cutoff"=>cutoff,"maxdim"=>maxdim)
end
