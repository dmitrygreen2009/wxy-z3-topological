# Infinite winding constraints: algebra passed, optimizer gate pending

The direct finite winding-QN implementation has passed its microscopic basis,
checkpoint-bridge and q0/q1 fixture gates for both cylinder families. Repeating
the same QN weights in an infinite MPS, by itself, fixes only a global product
charge and does not certify each microscopic winding loop.

An exact stronger construction is possible in principle. For zigzag, the
chosen winding operator has support inside one complete spatial slice. For
armchair it occupies two adjacent slices. Use nonoverlapping copies of these
operators, retaining their full matter action, and express their weights in
the validated local matter basis. Each exact Hamiltonian product term has
zero Z3 charge separately in every such winding cell it touches. This is
checked directly from the saved microscopic endpoint pairs for 20 available
star-ordered geometry/unit-cell combinations, widths 1–4 and cell sizes
2, 3 or 6 where recorded. No numerical eigensolver is used in that check.

`infinite_winding_charge_flow_audit.json` saves each explicit weight vector,
the exact source geometry, checked term counts, charge-period compatibility
and unavailable geometry records. The test uses the independently proved
full eight-state real exchange identity. Its two distinct number-flow types
are S_i^- σ^+ and S_i^+ S_j^- S_k^- σ^+; the controlled-Sz factor has the
same flow as the first type. Hermitian conjugates reverse zero flow and
therefore preserve the same charges. All terms preserve physical total N_up.

To constrain an individual winding-cell charge in an infinite MPS, virtual
charge support at each winding-cell boundary must be fixed explicitly, or
shown to remain restricted under every official solver operation. A recorded
constant background charge can make the chosen loop eigenvalue compatible
with zero virtual boundary flux. Physical number backgrounds must be treated
separately, with their actual rational filling reported. In a two-slice
zigzag cell, the two individual one-slice loops require separate accounting;
their product is insufficient. Armchair odd-slice source cells do not align
with the two-slice charge period and require an enlarged cell for this route.

No infinite constrained optimizer is certified yet. Initialization, subspace
expansion, the official VUMPS update and canonicalization must be checked for
boundary-block retention, then translated loop means and variances measured
on both geometries. Original infinite checkpoints remain unchanged. Existing
unrestricted energy/entropy jobs continue; these algebraic checks do not
invalidate them or restart them.

Fixing a disjoint loop-charge pattern may also constrain contractible CGS
combinations. It must not be called a single emergent flux sector or a MES
without checking those relations. This prerequisite establishes no
thermodynamic filling, converged entropy intercept, deconfinement, or 2D phase.

Physical individual matter occupations after the basis rotation are not mode
occupations. A future infinite measurement driver must either transform the
appropriate local matter operators or explicitly report only invariant star
number and physical gauge occupations. Full-state spatial entropy at cuts
between complete stars and complete-cell transfer eigenvalues are invariant
under the local matter unitary; that does not remove their convergence gates.
