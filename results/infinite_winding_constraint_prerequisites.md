# Infinite winding constraints: algebra and optimizer integration validated

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

A gated implementation now exists in `src/infinite_winding_basis.jl`. It retains all individual spin sites and all matter states, uses the exact real Hamiltonian, and initializes explicit physical-number backgrounds plus separate winding QNs for the disjoint loops in the cell. Winding background offsets are recorded separately from physical phase weights. The pinned library supports four QN fields, so this prototype permits up to three independent loop fields plus physical number. Each winding boundary is checked for zero charge support in both AL and AR.

`test/infinite_winding_basis.jl` tests both families and q0/q1 at a deliberately chosen fixture density 1/3, compares periodic microscopic energies in the original physical basis, verifies the full matter-dressed loop action, then checks both ordinary and translated loop values after canonicalization, expansion and one official infinite update. The test uses at most 18 individual sites, bond cap 3 and one update; it is not a cylinder convergence run. No production constrained infinite state is claimed before it passes. Existing checkpoints and calculations remain unchanged.

The initial gate passed all 110 checks in GitHub CI at ee75cdca: both geometries and q0/q1 retain their individual loop eigenvalues through canonicalization, expansion and an official infinite update. The extended gate passed 158 checks at 5112ec1, including the physical-species occupation identity. The production measurement/checkpoint wrapper passed seven integration checks in the same run. Downloaded machine-readable fixtures are saved in `infinite_winding_basis_fixture.json` and `infinite_winding_driver_ci_fixture.json`; selected log summaries are in `infinite_winding_driver_ci_test_summary.log`. Their density 1/3 and short updates are implementation fixtures, not optimized physical filling or phase evidence. The q2 fixture extension is a subsequent test change, not part of the 5112ec1 validation record.

The wrapper in `scripts/infinite_winding_qn.jl` requires an explicit physical N_up background and checks winding support and translated loop values after every official update. It preserves a QN optimization checkpoint for resume, separately from any dense measurement copy. Physical species magnetizations are reconstructed through the exact cyclic-symmetry identity; rotated mode occupations remain separately named. Source payload SHA is checked after measurement. No result is admitted to TEE fits by this driver.

A detected fixed-point rank of one from scalar Krylov starts is explicitly insufficient for automatic recanonicalization. The new measurement path requires a complete dense virtual spectrum (resource limit dimension 256, tolerance 1e-10) with one unit-modulus peripheral eigenvalue. At larger dimensions, or with an unresolved peripheral subspace, inconsistent centers are retained and their entropy is marked unresolved. Previously saved recanonicalized variational states remain preserved; any assertion of equivalence to an original nonprimitive boundary choice needs a separate complete source-spectrum audit. Their measured properties do not establish the physical phase.
