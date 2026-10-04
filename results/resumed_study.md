# Resumed study

The user accepted the audited full N_up=9 torus spectrum as the corrected
benchmark. Regressions now test the first ten ordered levels and multiplicities
4, 1, 1, 4. The independently audited next four levels remain regression checks.
`torus_hard_gate_audit.md` is preserved unchanged as the historical validation
record, including its then-current stop status.

Finite calculations resume from the saved individual-spin MPS checkpoints.
Star ordering has exactly the same spatial cut as the original axial ordering.
The renewed run uses single-threaded BLAS and no block-sparse multithreading;
the original multi-threaded attempt also enabled Strided permutation threads,
which the library warned could impair performance. Subsequent optional block
threading explicitly disables those extra Strided threads.

Infinite calculations now save the canonical state and measurements after
each dimension stage, rather than waiting until the entire calculation ends.
Transfer eigenvalues, residuals, canonical consistency, energy, full-state
spatial entropy, and magnetization are recorded. Values from an unconverged
state are preliminary; a finite transfer correlation length at one width or
bond dimension does not establish a finite ground-state correlation length.

Independent ED of the L=2, w=1 open cylinders checks the cylinder construction
in addition to the periodic torus. Zigzag uses 19 physical spins and a
92,378-dimensional N_up=10 sector; armchair uses 20 spins and a
184,756-dimensional N_up=10 sector. The independently obtained energies are
-7.696107978815031 and -7.994302046793368. ITensor errors are 5.01e-10 at
chi=128 and 1.86e-9 at chi=256, respectively. Energy variances are 3.55e-9 and
1.48e-8. The first chi=128 armchair comparison failed, triggering an increase
to 256; the validation thresholds were preserved. The zigzag full-state
spatial entropy also matches independent ED to 1.2e-8 nats. Armchair ED ground
states are degenerate, so a single ED eigenvector entropy is not a unique
reference for a DMRG-selected state.

The infinite-state canonical audit also measures left and right isometry
errors and center normalization. At zigzag w=1, chi=32 the isometry errors
are below 3e-15, while AL*C versus C*AR consistency is 0.00197. The transverse
order parameter mean |<S+>| is 0.280. Thus its transfer spectrum is a valid
spectrum of the recorded variational MPS, but its short correlation length is
not a converged ground-state conclusion. Further dimension refinement and
symmetry checks are required.

The L=4, w=1 star-ordering zigzag run found E=-15.0328573108 and
S=1.8713195070 at chi=256, below the original axial result
E=-15.0098830532 with S=0.8506911145. Its entropy changes by about 0.000257
between 128 and 256, with a further 4e-7 noiseless refinement change. This is
an explicit example of state selection affecting apparent entropy fits, not
just small truncation corrections. Multiple seeds and additional lengths
are included in the continuation workflow. Width one also has parallel-edge
identifications; a fourth circumference is included so future fits can omit
that special narrow quotient.

The initial U(1)-conserving product-state VUMPS attempt encountered the
library's documented zero-norm two-site expansion for long-range QN
Hamiltonians. Star ordering improves the expansion, but a low-dimension
state can still have near-machine canonical residuals and a high energy with
frozen bonds. Those trial states are retained as algorithm diagnostics and
must not be used as ground-state correlation lengths. A revised product
initialization activates nearest matter-gauge hopping pairs at both A and B
stars. An independent alternative uses the official infinitemps_approx
variational fit of a finite DMRG state, followed by canonicalization and
VUMPS. Transfer matrices are never formed by arbitrary identification of
unrelated finite bond spaces.

Repository synchronization checkpoint: validated source, corrected benchmark
regressions, full geometry records, environments, pipeline validation, and the
first resumed finite grid were committed and pushed to `main` at
`5dfb1876638f8a38032c27f821174fc069849070`. This is a reproducibility milestone,
not a scientific endpoint. Legacy provenance limitations are explicit in
`legacy_provenance_audit.json`; recoverable iteration values are in
`recorded_iterations.csv`.

The self-contained audited VUMPS driver regression passes both tests without
requiring an untracked MPS. The active finite batch adds width four, then checks
higher bond dimensions and longer cylinders. Independent seed runs test the
large width-one entropy change accompanying a lower variational energy. No
thermodynamic intercept is claimed from the current short, unconverged grid.

Independent seed audit on the star-ordered zigzag L4 width-one cylinder at
chi=256: seeds 7103 and 7104 reproduce E=-15.0328573108151 and
S=1.871319507. All four parallel-edge two-cycle expectations are nearly one.
Seed 7105 instead reaches E=-15.0098845261812, S=0.850744667, with nontrivial
loop content in the last two cycles. Direct variances are 6.53e-8 and 1.04e-8
for seeds 7104 and 7105 respectively. The latter is a sector-selection failure
for a ground-state entropy fit despite its small sweep drift and variance.
Raw points and loop/variance diagnostics are retained in `seed_selection_audit.json`.
This is finite microscopic evidence about solver state selection, not evidence
establishing topological order.

Further validation checkpoint: independent ED across every number sector confirms
the global ground sectors [3], [5,6], and [9] for the star, pair, and torus.
This scan does not estimate degeneracies; the separate block spectrum audit
remains the multiplicity validation. Exact one-particle spectra on nine explicit
geometries agree with the singular values implied by W†W=I to 3.22e-15.
Neither observation establishes the interacting thermodynamic phase.

The independent infinite fixed-point canonicalization pipeline passes analytic
and QN entangled-state tests, including originally untagged, primed virtual
indices. Measurement canonicalization preserves the optimized state and its
solver residual. Zigzag width-one chi32 and chi64 yield full-state entropies
0.9664895127 and 1.1224124645, and correlation lengths in axial slices
2.7855843155 and 6.8977827251. The strong bond-dimension dependence precludes
a converged finite-correlation-length claim. Leading transfer eigenvalues and
canonicalization errors are saved in the fixedpoint JSON records.

Exact local CGS projection supplies a better seed without replacing the
microscopic Hamiltonian. Projected zigzag width-one chi512 results at L4, L6,
and L10 have entropies 1.8713200977, 1.7446134177, and 2.1079932659. L4
agrees with chi256 within 5.4e-7; L6 still has appreciable energy drift and
L10 remains limited by truncation. These lengths cannot yet establish an
asymptotic entropy. Both width-four L4 chi256 points completed, but truncation
and refinement changes make them exploratory fit inputs only. Raw data and
fits remain separate; fitted intercepts are not topological conclusions.

Independent eight-vector block ED on the open L2 width-one cylinders finds
one ground state for zigzag (fixed-sector gap 0.04434958718) and four for
armchair (gap 0.01252816539). Maximum residuals are 2.61e-13 and 5.24e-13.
This again shows why a scalar low-eigenvalue list cannot certify degeneracies:
the latest armchair scalar run returned only two of the four ground states.
These tiny open-cylinder gaps are not thermodynamic gap estimates.
