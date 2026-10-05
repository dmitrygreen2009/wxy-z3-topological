# Sector-preserving optimization at independently audited finite fillings

Checkpoint classification: VALIDATED AND USABLE FOR PHASE INFERENCE: none
yet. COMPUTED BUT NOT YET CONVERGED: the wider-cylinder continuations.
REJECTED / DIAGNOSTIC ONLY: the accepted small zigzag and armchair q0/q1/q2
optimizer benchmarks and exact representation/initializer
fixtures (usable for implementation validation, awaiting 2D gates), plus the
preserved failed cold initialization. The q0 state is retained; this category
does not reject its numerical values.

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

The armchair counterparts are now complete as recorded below. The already validated N_up=2 fixture comparisons establish
sector-algorithm behavior, not physical filling. No point here is admitted to
a topological-entropy fit: common MES identification and length/circumference
convergence remain separate requirements. No finite-ring substitution for an
infinite transfer matrix is used, and no finite-state correlation length is
assigned by joining unrelated end bonds. No 2D phase is established or excluded.

## Completed zigzag physical-filling winding-sector candidates

All three separately optimized q=0,1,2 states pass the frozen fixed-number /
winding-QN gates. They share N_up=9, S^z=-1/2, filling 9/19, 19 physical
spins and the full-state spatial cut at bond 10. Number and winding charge
are imposed exactly by the MPS QNs; they are not chosen variationally.
The full finite-geometry number scan supplies the filling comparison.

| q | Energy | ED error | H variance | Full-state S | Cap / achieved bond | Loop purity error |
|---|---:|---:|---:|---:|---:|---:|
| 0 | -7.696107978803282 | 1.177e-11 | 8.372e-11 | 1.355308930211669 | 107 / 107 | 1.074e-15 |
| 1 | -7.651758391620112 | 1.153e-11 | 7.159e-11 | 0.867977280788276 | 128 / 113 | 8.327e-16 |
| 2 | -7.651758391620521 | 1.111e-11 | 6.871e-11 | 0.867977280823546 | 128 / 113 | 3.511e-16 |

The q1/q2 runs progressed through caps 32,64,128. Cap64 failed the
unchanged energy-variance and ED-error gates despite small sweep drift;
cap128 passed, with nine final zero-noise sweeps and entropy changes
9.79e-8 / 9.64e-8 from the cap64 states. The q1/q2 winding expectations
are omega / conjugate(omega), respectively, with loop variances consistent
with zero at roundoff. Their energy splittings above q0 are
0.0443495871831705 and 0.0443495871827608 J. This width-one splitting is
not a circumference-scaling result or a confinement measurement.

The complete scalar-ARPACK sector lists are not multiplicity certificates.
No uniqueness or MES claim follows from these small excited-sector
solutions. The method, exact full eight-state matter representation,
Hamiltonian MPO and winding constraints are frozen after the independent
ED and regression validations; the later exact total-spin basis audit
identifies the same basis up to signed permutation. Older initializer
reachability restrictions do not invalidate these converged results, which
pass independent ED and variance checks. The q1/q2 run process retained
its launch implementation and did not reload later source edits.

Correlation lengths remain null for these finite candidates: no infinite
transfer matrix can be obtained by joining their distinct open end bonds.
They remain REJECTED / DIAGNOSTIC ONLY for 2D phase interpretation, while
valid and retained as quantitative sector-preserving optimizer benchmarks.
No entropy point has been admitted to a gamma fit.

## Completed armchair physical-filling winding-sector candidates

All three separately optimized armchair L2 width-one candidates pass the
frozen fixed-number/winding gates. Each retains all 20 physical spins with
fixed N_up=8, S^z=-2, filling 8/20, and the validated full-state spatial
cut recorded in its result JSON. The independent all-number finite scan
identifies N_up=8 and its particle-hole partner N_up=12 as global finite
minima; it does not determine the infinite filling.

| q | Energy | ED error | H variance | Full-state S | Cap / achieved bond | Loop purity error |
|---|---:|---:|---:|---:|---:|---:|
| 0 | -8.037741645043154 | 1.439e-11 | 1.311e-10 | 1.144144154081445 | 198 / 198 | 1.385e-15 |
| 1 | -7.495281096629017 | 3.524e-12 | 3.208e-11 | 1.247446159039853 | 256 / 224 | 1.535e-15 |
| 2 | -7.495281096627773 | 4.780e-12 | 4.493e-11 | 1.338533480820499 | 256 / 222 | 9.550e-16 |

The q0 source is the preserved physical-spin N_up8 ground checkpoint,
converted through the exact basis/projector bridge and refined in the
winding QN representation. Nontrivial q1/q2 branches were independently
optimized through caps 32,64,128,256; insufficient earlier caps remain
recorded. Their winding means are omega and conjugate(omega), and loop
variances vanish within roundoff. Their spectra contain independently
observed degenerate ground vectors within the fixed q sector. Different
entropies in q1 and q2 therefore do not establish different quantum
dimensions or a failure of energy convergence: these optimized branches
have not been certified as minimally entangled states.

All six small physical-filling sector candidates now validate the exact
sector-preserving implementation against independent ED. No complete
multiplicity claim is made from scalar ARPACK. Correlation lengths remain
null for these open finite states, and all six remain REJECTED / DIAGNOSTIC
ONLY for 2D phase inference pending the distinct MES, length, circumference
and bulk-filling gates. No point is admitted to an entropy-intercept fit.

At the user's resource checkpoint, 12 lower-priority Julia jobs were
confirmed against their purposes and disk checkpoints, then paused with
SIGSTOP. The active armchair sector batch was preserved and completed
normally. A lightweight watcher successfully resumed the first queued
width-two armchair number-sector refinement with SIGCONT; eleven jobs
remain paused and the active heavy-Julia limit is one. Checkpoints, parent
processes, and in-memory progress were preserved. The machine-local queue
and independent checksum/process verification are recorded separately in
`julia_resource_pause_audit.json` and `julia_resource_pause_verification.json`.
