# Rebuild the infinite state's Schmidt centers from its AL transfer fixed point
# using official library canonicalization, without any variational update.
function canonicalize_left(AL;tol=1e-14)
    # Densification preserves amplitudes and is confined to measurement.
    hasqns(siteind(AL,1)) && (AL=InfiniteMPS([dense(AL[j]) for j=1:nsites(AL)],translator(AL)))
    # Official canonicalization expects unprimed links explicitly tagged Left.
    # Some VUMPS snapshots use untagged links and prime levels instead.
    AL=InfiniteMPS([addtags(noprime(AL[j]),"Left";tags="Link") for j=1:nsites(AL)],translator(AL))
    _,right,lambda_right=ITensorInfiniteMPS.right_orthogonalize(AL;left_tags=ts"Left",right_tags=ts"Right",tol)
    left,centers,lambda_left=ITensorInfiniteMPS.left_orthogonalize(right;tol)
    @assert abs(lambda_left-1)<1e-10
    InfiniteCanonicalMPS(left,centers,right),lambda_right
end
function infinite_schmidt_entropies(psi)
    map(1:nsites(psi)) do j
        probabilities=abs2.(svdvals(Array(psi.C[j],inds(psi.C[j])...)))
        probabilities./=sum(probabilities)
        -sum(p>0 ? p*log(p) : 0.0 for p in probabilities)
    end
end
