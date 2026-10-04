# Local microscopic symmetry sectors in infinite width-one states

Independent microscopic loop-MPO tests pass: U^3=I on an arbitrary complex
state, the exactly known product-state charge, and the charge-three matter
operator expectation. The canonicalized QN zigzag width-one chi32 and chi64
states have Re<U> approximately -1/2 at each of the two cell slices, but
Im<U> is far from either ±sqrt(3)/2. Their local cycle charges are mixed.
These values are saved in the fixedpoint observable JSON files.

The tiny zigzag quotient has an additional **local S3 symmetry**. Exchange
its two parallel physical gauge edges (legs 1 and 3 at both endpoints),
interchange matter species a with -a mod 3, and multiply each matter ket by
omega^(-sum_a a*n_up[a]), with a=0,1,2. This reflection R commutes with both
endpoint Hamiltonians, including their unchanged external-leg couplings.
It obeys R²=I, U³=I, and RUR=U†. Independent ED verifies every identity in
all eleven number blocks of the complete 1024-ket two-star, ten-spin space;
see `parallel_edge_symmetry_audit.json` and the explicit L1 width-one graph.
Thus nontrivial local q=1 and q=2 charges form degenerate doublets whenever
that local representation occurs. This statement does not assume that the
global ground state belongs to that representation. Wider zigzag quotients
lack this particular parallel-edge exchange symmetry.

Consequently, the measured width-one entropy and transfer-length growth
cannot automatically be interpreted as a bulk critical mode or a
thermodynamic topological constant. Mixed local representations and state
selection require investigation. The finite projected L4 and refined L6
width-one states instead have every local two-cycle expectation nearly one,
as independently measured with the full microscopic gate action. Their
entropies must not be merged with the infinite mixed-cycle states as though
they were identical sectors.

The charge-three matter one-point expectation vanishes in these QN states,
as required by conserved U1; this does not rule out spontaneous matter
order. Full-state entropies still include every matter and gauge spin.
No gauge-only entropy is used as topological evidence.
