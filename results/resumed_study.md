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
