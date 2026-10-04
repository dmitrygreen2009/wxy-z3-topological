# Microscopic observables evaluated by the library's InfiniteCanonicalMPS expect.
isdefined(@__MODULE__, :cgs_product_mpo) || include("cgs_operator_mpo.jl")
function infinite_twocycle_mpo(s,offset=0)
    # Scalar translated indices remain valid beyond a CelledVector unit cell.
    cgs_product_mpo(s,cgs_circumference_spec("zigzag",1;offset))
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
