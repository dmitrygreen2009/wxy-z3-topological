# Domain-specific initialization only. A deliberate unit-dimensional cell
# boundary repeats independent finite-cell states; no unrelated entangled
# finite bonds are identified. Production remains on individual spin sites.
function infinite_cell_product(phi::MPS;tol=1e-14)
    @assert !hasqns(phi) || iszero(flux(phi)) "Use the library's shifted cell charge convention"
    @assert abs(norm(phi)-1)<1e-12
    n=length(phi)
    boundary=hasqns(phi) ? Index([zero(flux(phi))=>1],"Link,Left,l=$n";dir=ITensors.Out) : Index(1,"Link,Left,l=$n")
    boundary=addtags(boundary,ITensorInfiniteMPS.celltags(1))
    previous=ITensorInfiniteMPS.translatecelltags(boundary,-1)
    nodes=[addtags(addtags(phi[j],"Left";tags="Link"),ITensorInfiniteMPS.celltags(1)) for j=1:n]
    nodes[1]*=onehot(previous=>1);nodes[end]*=onehot(dag(boundary)=>1)
    A=InfiniteMPS(nodes,ITensorInfiniteMPS.translatecelltags)
    @assert dim(only(linkinds(A,n=>n+1)))==1
    _,right,_=ITensorInfiniteMPS.right_orthogonalize(A;left_tags=ts"Left",right_tags=ts"Right",tol)
    left,center,lambda=ITensorInfiniteMPS.left_orthogonalize(right;tol)
    @assert abs(lambda-1)<1e-10
    # Official polar canonicalization retains the input arrow at the special
    # one-dimensional boundary. Right-canonical subspace expansion expects the
    # opposite arrow. Reorient this zero-QN index consistently in AR and C;
    # storage, amplitudes, QN labels and the physical state are unchanged.
    if hasqns(phi)
        i=only(commoninds(right[n],right[n+1]))
        if ITensors.dir(i)!=ITensors.Out
            @assert dim(i)==1 && iszero(qn(i,ITensors.Block(1)))
            flip(T,k)=ITensors.setinds(T,[j==k ? dag(j) : j for j in inds(T)])
            right[n]=flip(right[n],i)
            right[n+1]=flip(right[n+1],i)
            center[n]=flip(center[n],i)
        end
    end
    InfiniteCanonicalMPS(left,center,right)
end
