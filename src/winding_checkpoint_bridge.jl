# Domain-specific exact operator only. ITensorMPS constructs/decomposes/applies
# every MPO; this file does not implement a tensor-network solver/contraction.
isdefined(@__MODULE__,:matter_charge_basis) || include("matter_charge_basis.jl")
function inverse_matter_basis_mpo(lat,s;cutoff=1e-14)
    result=MPO(ComplexF64,s);bonds=[linkind(result,j) for j=1:lat.n-1]
    for j=1:lat.n
        T=op("Id",s[j]);j>1 && (T*=onehot(dag(bonds[j-1])=>1));j<lat.n && (T*=onehot(bonds[j]=>1));result[j]=T
    end
    ip=invperm(lat.order);B=matter_charge_basis()[8:-1:1,8:-1:1]'
    for v=1:lat.nv
        positions=[ip[3(v-1)+a] for a=1:3]
        @assert positions==collect(first(positions):first(positions)+2)
        ss=s[positions];gate=ITensor(reshape(B,2,2,2,2,2,2),prime.(ss)...,dag.(ss)...)
        block=MPO(gate,ss;cutoff);firstsite=first(positions)
        for k=1:3
            T=block[k]
            k==1 && firstsite>1 && (T*=onehot(dag(bonds[firstsite-1])=>1))
            k==3 && firstsite+2<lat.n && (T*=onehot(bonds[firstsite+2]=>1))
            result[firstsite+k-1]=T
        end
    end
    result
end
function lifted_winding_link(old,j)
    # Library SVD links can point In. Preserve their signed number flow when
    # expressing every lifted link with the Out orientation.
    sign=hasqns(old) && dir(old)==ITensors.In ? -1 : 1
    blocks=hasqns(old) ? [(sign*val(qn(old,ITensors.Block(b)),"Sz"),blockdim(old,ITensors.Block(b))) for b=1:nblocks(old)] : [(0,dim(old))]
    !hasqns(old) && @assert dim(old)==1 "A nontrivial dense operator link lacks its number grading"
    space=Pair{QN,Int}[];mapping=zeros(Int,dim(old),3);oldoffset=0;offset=0
    for (sz,d) in blocks
        for q=0:2
            push!(space,QN(("Sz",sz),("Winding",q,3))=>d)
            for k=1:d;mapping[oldoffset+k,q+1]=offset+k;end
            offset+=d
        end
        oldoffset+=d
    end
    Index(space,"Link,WindingBridge,l=$j";dir=ITensors.Out),mapping
end
function winding_basis_projector_mpo(lat,oldsites,p;charge,cutoff=1e-14)
    @assert charge in 0:2
    @assert all(hasqns(s) && val(qn(s,ITensors.Block(1)),"Sz")==1 && val(qn(s,ITensors.Block(2)),"Sz")==-1 for s in oldsites)
    weights=winding_site_weights(lat,p);newsites=winding_qn_sites(weights)
    R=inverse_matter_basis_mpo(lat,oldsites;cutoff)
    oldlinks=[linkind(R,j) for j=1:lat.n-1]
    lifted=[lifted_winding_link(i,j) for (j,i) in enumerate(oldlinks)]
    newlinks=[x[1] for x in lifted];maps=[x[2] for x in lifted];nodes=ITensor[]
    for j=1:lat.n
        originalinds=Index[prime(oldsites[j]),dag(oldsites[j])]
        j>1 && push!(originalinds,dag(oldlinks[j-1]));j<lat.n && push!(originalinds,oldlinks[j])
        dL=j==1 ? 1 : dim(oldlinks[j-1]);dR=j==lat.n ? 1 : dim(oldlinks[j])
        values=reshape(Array(R[j],originalinds...),2,2,dL,dR)
        newdL=j==1 ? 1 : dim(newlinks[j-1]);newdR=j==lat.n ? 1 : dim(newlinks[j])
        liftedvalues=zeros(ComplexF64,2,2,newdL,newdR)
        for out=1:2,input=1:2,left=1:dL,right=1:dR
            value=values[out,input,left,right];iszero(value) && continue
            for qleft in (j==1 ? (0:0) : (0:2))
                qright=mod(qleft-(out==1 ? weights[j] : 0),3)
                j==lat.n && qright!=mod(-charge,3) && continue
                l=j==1 ? 1 : maps[j-1][left,qleft+1]
                r=j==lat.n ? 1 : maps[j][right,qright+1]
                liftedvalues[out,input,l,r]=value
            end
        end
        newinds=Index[prime(newsites[j]),dag(oldsites[j])]
        j>1 && push!(newinds,dag(newlinks[j-1]));j<lat.n && push!(newinds,newlinks[j])
        push!(nodes,ITensor(reshape(liftedvalues,dim.(newinds)...),newinds...))
    end
    operator=MPO(nodes)
    @assert val(flux(operator),"Sz")==0 && mod(val(flux(operator),"Winding"),3)==charge
    (;operator,newsites,weights,charge,interpretation="Exact P_q B†: all original matter states retained, onsite winding projected, number preserved")
end
function bridge_physical_checkpoint(psi,lat,p;charge,cutoff=1e-13,maxdim=4maxlinkdim(psi))
    @assert hasqns(psi)
    oldnumber=val(flux(psi),"Sz")
    bridge=winding_basis_projector_mpo(lat,siteinds(psi),p;charge,cutoff=1e-14)
    converted=apply(bridge.operator,psi;cutoff,maxdim)
    weight=norm(converted)^2/norm(psi)^2
    @assert weight>1e-14 "Requested winding sector has negligible weight; preserve the source and choose another initializer"
    normalize!(converted)
    @assert val(flux(converted),"Sz")==oldnumber && mod(val(flux(converted),"Winding"),3)==charge
    converted,Dict("projection_weight"=>weight,"charge"=>charge,"number_preserved"=>true,"operator_maxlinkdim"=>maxlinkdim(bridge.operator),
        "conversion_cutoff"=>cutoff,"conversion_maxdim"=>maxdim,"basis"=>"exact_matter_charge_basis","onsite_winding_weights"=>bridge.weights)
end
