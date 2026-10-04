# Microscopic observables evaluated by the library's InfiniteCanonicalMPS expect.
function infinite_twocycle_mpo(s,offset=0)
    @assert length(s)>=offset+8
    u=zeros(ComplexF64,8,8)
    # p=(1,0,2) at both endpoints: c=1 and matter permutation q=2.
    for old=0:7
        new=sum(((old>>a)&1)<<mod(a+2,3) for a=0:2)
        u[new+1,old+1]=cis(2pi*(3-count_ones(old))/3)
    end
    gate(js)=ITensor(reshape(u,2,2,2,2,2,2),prime.(s[js])...,dag.(s[js])...)
    U=gate(offset.+[1,2,3])*gate(offset.+[6,7,8])
    for (j,p) in [(offset+4,1),(offset+5,2)]
        U*=ITensor(ComplexF64[cis(2pi*p/3) 0;0 1],s[j]',dag(s[j]))
    end
    MPO(U,collect(s[offset+1:offset+8]);cutoff=1e-14)
end
function infinite_triplet_mpo(s,first)
    os=OpSum();os+=1.0,"S+",1,"S+",2,"S+",3
    MPO(os,collect(s[first:first+2]))
end

# CGS invariant, U(1) charge-three matter correlation. Both triplets retain
# their physical matter spins; intervening sites carry the identity.
function infinite_triplet_pair_mpo(s,first,second)
    @assert second>=first+3
    sites=[s[j] for j=first:second+2]
    os=OpSum()
    r=second-first+1
    os+=1.0,"S+",1,"S+",2,"S+",3,"S-",r,"S-",r+1,"S-",r+2
    MPO(os,sites)
end
