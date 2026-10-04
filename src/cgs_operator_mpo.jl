# Domain-specific microscopic observable construction. All decomposition,
# application, and contraction use ITensor/ITensorMPS library routines.
# Factor into local three-spin permutation MPOs; never form a circumference
# operator as a dense tensor with exponentially many matrix entries.
function cgs_circumference_spec(family,w;offset=0)
    @assert family in ["zigzag","armchair"] && w>=1
    triplets=Tuple{Int,Int,Int}[];gauges=Tuple{Int,Int}[]
    for j=0:w-1
        row=offset+9j
        if family=="zigzag"
            append!(triplets,[(row+1,1,2),(row+6,1,2)])
            append!(gauges,[(row+4,1),(row+5,2)])
        else
            append!(triplets,[(row+1,1,2),(row+6,1,1),(row+9w+1,1,1),(row+9w+6,1,2)])
            append!(gauges,[(row+4,1),(row+5,2),(row+9,2),(row+9w+4,1)])
        end
    end
    first_site=minimum(first.(triplets));last_site=maximum(vcat([j+2 for (j,c,q) in triplets],first.(gauges)))
    (;family,w,offset,triplets,gauges,first_site,last_site)
end

function cgs_disjoint_pair_spec(first,second)
    # The second loop is adjointed. These disjoint-support products use the
    # exact Abelian microscopic representation U(c,q)^dagger=U(-c,-q).
    support_first=Set(vcat([j+k for (j,c,q) in first.triplets for k=0:2],[j for (j,p) in first.gauges]))
    support_second=Set(vcat([j+k for (j,c,q) in second.triplets for k=0:2],[j for (j,p) in second.gauges]))
    @assert isempty(intersect(support_first,support_second))
    triplets=vcat(first.triplets,[(j,mod(-c,3),mod(-q,3)) for (j,c,q) in second.triplets])
    gauges=vcat(first.gauges,[(j,mod(-p,3)) for (j,p) in second.gauges])
    (;triplets,gauges,first_site=min(first.first_site,second.first_site),last_site=max(first.last_site,second.last_site))
end
function cgs_product_mpo(s,spec;cutoff=1e-14)
    # CelledVector translates scalar indices beyond its stored unit cell, but
    # generic range indexing applies finite-array bounds checks.
    ss=[s[j] for j=spec.first_site:spec.last_site];N=length(ss)
    result=MPO(ComplexF64,ss);bonds=[linkind(result,j) for j=1:N-1]
    @assert all(dim(b)==1 for b in bonds)
    # The product-operator constructor truncates/canonicalizes its identity
    # tensors. Those tensors can carry compensating scalar phases; replacing
    # selected tensors would lose those phases. Fill an explicit local-identity
    # scaffold instead, using the library's unit-bond indices and operators.
    for j=1:N
        node=op("Id",ss[j])
        j>1 && (node*=onehot(dag(bonds[j-1])=>1))
        j<N && (node*=onehot(bonds[j]=>1))
        result[j]=node
    end
    for (first,c,q) in spec.triplets
        u=zeros(ComplexF64,8,8)
        for old=0:7
            new=sum(((old>>a)&1)<<mod(a+q,3) for a=0:2)
            u[new+1,old+1]=cis(2pi*c*(3-count_ones(old))/3)
        end
        local_sites=[s[j] for j=first:first+2]
        tensor=ITensor(reshape(u,2,2,2,2,2,2),prime.(local_sites)...,dag.(local_sites)...)
        block=MPO(tensor,local_sites;cutoff)
        rel=first-spec.first_site+1
        for k=1:3
            node=block[k]
            k==1 && rel>1 && (node*=onehot(dag(bonds[rel-1])=>1))
            k==3 && rel+2<N && (node*=onehot(bonds[rel+2]=>1))
            result[rel+k-1]=node
        end
    end
    for (j,p) in spec.gauges
        rel=j-spec.first_site+1;node=ITensor(ComplexF64[cis(2pi*p/3) 0;0 1],s[j]',dag(s[j]))
        rel>1 && (node*=onehot(dag(bonds[rel-1])=>1))
        rel<N && (node*=onehot(bonds[rel]=>1))
        result[rel]=node
    end
    result
end
