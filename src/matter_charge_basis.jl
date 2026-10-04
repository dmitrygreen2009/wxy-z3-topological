# Exact eight-state local basis rotation, not a matter/clock truncation.
function matter_charge_basis()
    B=zeros(ComplexF64,8,8);B[1,1]=B[8,8]=1
    for a=0:2,b=0:2
        B[(1<<b)+1,(1<<a)+1]=cis(-2pi*a*b/3)/sqrt(3)
        B[(7 ⊻ (1<<b))+1,(7 ⊻ (1<<a))+1]=cis(2pi*a*b/3)/sqrt(3)
    end
    @assert norm(B'B-I)<1e-13
    B
end
function winding_site_weights(lat,p)
    @assert length(p)==lat.n-3lat.nv
    @assert all(mod(sum(p[e] for e in es),3)==0 for es in lat.legs)
    weights=zeros(Int,lat.n);ip=invperm(lat.order)
    for v=1:lat.nv
        c=mod(p[lat.legs[v][1]],3);q=mod(p[lat.legs[v][2]]-c,3)
        for a=0:2;weights[ip[3(v-1)+a+1]]=mod(c+q*a,3);end
    end
    for e in eachindex(p);weights[ip[3lat.nv+e]]=mod(p[e],3);end
    weights
end
function winding_qn_sites(weights)
    [Index([QN(("Sz",1),("Winding",w,3))=>1,QN(("Sz",-1),("Winding",0,3))=>1],"Site,S=1/2,n=$j") for (j,w) in enumerate(weights)]
end
function winding_initial_state(weights,nup,charge)
    n=length(weights);@assert 0<=nup<=n && charge in 0:2
    feasible=falses(n+1,nup+1,3);feasible[1,1,1]=true
    for j=1:n,k=0:min(j,nup),q=0:2
        feasible[j+1,k+1,q+1]=feasible[j,k+1,q+1] || (k>0 && feasible[j,k,mod(q-weights[j],3)+1])
    end
    @assert feasible[n+1,nup+1,charge+1] "Requested number/winding block is empty"
    state=fill("Dn",n);k=nup;q=charge
    for j=n:-1:1
        down=feasible[j,k+1,q+1]
        up=k>0 && feasible[j,k,mod(q-weights[j],3)+1]
        if up && (!down || rand(Bool))
            state[j]="Up";k-=1;q=mod(q-weights[j],3)
        end
    end
    @assert k==0 && q==0
    state
end
function charge_basis_opsum(lat;roundoff_cutoff=1e-14)
    B=matter_charge_basis();os=OpSum();ip=invperm(lat.order)
    for i=1:3
        A=zeros(ComplexF64,8,8)
        for a=0:2,old=0:7
            ((old>>a)&1)==1 || continue
            A[(old ⊻ (1<<a))+1,old+1]+=W[a+1,i]
        end
        rotated=B'A*B
        for old=0:7,out=0:7
            coefficient=rotated[out+1,old+1]
            abs(coefficient)>roundoff_cutoff || continue
            @assert count_ones(out)==count_ones(old)-1
            Q(b)=sum(a*((b>>a)&1) for a=0:2)
            @assert mod(Q(out)-Q(old)+i-1,3)==0
            for v=1:lat.nv
                term=Any[-coefficient];hc=Any[-conj(coefficient)]
                for a=0:2
                    input=(old>>a)&1;output=(out>>a)&1
                    opname=input==output ? (input==1 ? "ProjUp" : "ProjDn") : (output==1 ? "S+" : "S-")
                    hcname=input==output ? opname : (input==1 ? "S+" : "S-")
                    j=ip[3(v-1)+a+1];append!(term,[opname,j]);append!(hc,[hcname,j])
                end
                j=ip[3lat.nv+lat.legs[v][i]]
                append!(term,["S+",j]);append!(hc,["S-",j]);add!(os,term...);add!(os,hc...)
            end
        end
    end
    os
end
function physical_from_charge_basis(psi,lat;cutoff=1e-14,maxdim=4maxlinkdim(psi))
    physical=dense(psi);sites=siteinds(physical);ip=invperm(lat.order)
    B=matter_charge_basis()[8:-1:1,8:-1:1];gates=ITensor[]
    for v=1:lat.nv
        ss=[sites[ip[3(v-1)+a]] for a=1:3]
        push!(gates,ITensor(reshape(B,2,2,2,2,2,2),prime.(ss)...,dag.(ss)...))
    end
    apply(gates,physical;cutoff,maxdim)
end
