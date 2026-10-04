# Exact antiunitary particle-hole map Θ = (product of spin flips) K.
# K conjugates W and spin flips interchange each hopping term with its adjoint.
# Dagger every QN tensor, then pair its incoming physical index with a second
# outgoing index through the zero-flux spin-flip tensor. This preserves QN MPS
# structure while reversing total Sz, using only library tensor operations.
function particle_hole_mps(psi)
    s=siteinds(only,psi)
    tensors=[noprime(ITensor(ComplexF64[0 1;1 0],s[j],s[j]')*dag(psi[j])) for j=1:length(psi)]
    MPS(tensors)
end
