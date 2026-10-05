# Sector-preserving optimization at independently audited finite fillings

The zigzag L2 width-one trivial winding state is now converged under the
unchanged direct-QN acceptance criteria. It retains all 19 physical spins,
including local matter, in the exact microscopic Hamiltonian. Its fixed
physical number is N_up=9, total S^z=-1/2, filling 9/19. The complete preserved
small-cylinder number scan independently identifies N_up=9 and 10 as global
minima for this finite geometry. This does not select an infinite-system
filling.

| Quantity | Value |
|---|---:|
| Exact winding eigenvalue | 1 (q=0) |
| MPS energy | -7.696107978803282 |
| Independent fixed-sector ED energy | -7.696107978815055 |
| Energy variance | 8.3716145127255e-11 |
| Full-state spatial entropy | 1.355308930211669 |
| Entropy change between zero-noise passes | 7.549516567451064e-15 |
| Loop purity error | 1.0739257495113329e-15 |
| Loop variance | 0 |
| Bond cap and achieved dimension | 107 |

The spatial cut is bond 10 of the saved star ordering, with complete matter
triplets on each side. The exact local matter basis rotation therefore
preserves this entropy. The run uses 24 sweeps in two zero-noise passes,
cutoff 1e-13, Krylov dimension 12, eigensolver tolerance 1e-13, maximum 30
eigensolver iterations, and seed 7254. Its launch code is 5112ec1; exact
executed driver/dependency hashes, runtime, every sweep, Schmidt probabilities,
geometry references and immutable checkpoint SHA are recorded in
`zigzag_L2_w1_Nup9_winding0_qn_seed7254.json`.

Initialization uses the previously validated physical spin checkpoint,
unchanged on disk. The exact antiunitary particle-hole symmetry maps its
N_up=10 state to N_up=9. The exact P_0 B† bridge then maps it to the
number/winding QN basis. Projection weight is 0.9999999999906488; energy and
entropy preservation checks pass without loosening tolerances. Conversion
provenance and source SHA are in
`zigzag_L2_w1_Nup9_winding0_basis_bridge.json`.

The earlier cold initialization is preserved in
`archive/winding_optimizer_plateau/zigzag_L2_w1_Nup9_winding0_cold_plateau.json`.
It is rejected: tiny truncation error and sweep drift coexisted with variance
3.36e-5 and an ED energy error 6.44e-6. The successful alternative seed does
not establish the precise mechanism of that optimizer trap.

Nontrivial physical-filling winding states and the armchair counterparts
remain in progress. The already validated N_up=2 fixture comparisons establish
sector-algorithm behavior, not physical filling. No point here is admitted to
a topological-entropy fit: common MES identification and length/circumference
convergence remain separate requirements. No finite-ring substitution for an
infinite transfer matrix is used, and no finite-state correlation length is
assigned by joining unrelated end bonds. No 2D phase is established or excluded.
