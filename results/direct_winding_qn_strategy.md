# Direct microscopic winding quantum numbers

An exact local basis change diagonalizes all nine affine CGS matter actions
simultaneously. Let B be identity on the zero- and three-up matter states.
For one up spin at a, use
B|a> = (1/sqrt(3)) sum_b omega^(-ab)|b>.
For two up spins with the hole at a, use
B|hole a> = (1/sqrt(3)) sum_b omega^(ab)|hole b>.
Indices a,b are 0,1,2. Every star retains all eight microscopic matter states.
B commutes with total matter N_up, and is unitary. This is an exact spin-basis
change, with no clock approximation or low-energy projection.

For a vertex action U(c,q),
B† U(c,q) B = omega^[sum_a (c+qa)n_a].
The gauge factor remains omega^(p_e n_e). Thus a specified full microscopic
CGS winding operator becomes a product of onsite phases, with each matter
binary site's Z3 weight c+qa and each physical gauge site's weight p_e.
Each site still has dimension two. ITensor's U1 x Z3 block structure can
therefore fix total physical N_up and winding charge together. Local basis
weights depend on the chosen explicit cycle, not an assumed effective Wilson
loop. The relation of CGS eigenvalues to emergent topological flux remains a
separate physical identification.

The exact transformed exchange operator for leg i is
B† (sum_a W_ai mu^-_a) B. It has five nonzero matrix units in the full eight-state
matter space for each leg. Each lowers matter number by one and changes its
mode charge by -(i-1) mod 3. Multiplying by the incident physical sigma^+
conserves both onsite charge and total physical N_up. The MPO is built using
library OpSum/MPO from three individual matter sites and one shared gauge
site per term, with Hermitian conjugates. No whole-circumference supersite or
custom tensor-network contraction is introduced.

Independent audits:

- `results/matter_charge_basis_audit.json`: all eight states, unitarity,
  number preservation, all nine CGS actions, and exact exchange selection rules.
- `results/charge_basis_ed_audit.json`: full N_up=2 microscopic blocks on both
  smallest test cylinders, exact matrix conjugation of H, direct winding block
  selection, and agreement with original-basis sector ED. New block spectra
  and multiplicities are saved. This validates the representation, not filling.
- `test/matter_charge_basis.jl`: production-library QNs, original physical-H
  energies, microscopic loop eigenvalues, physical filling, and invariance
  of the full-state spatial Schmidt entropy. All 42 checks passed.

The unitary acts within each complete star. All audited spatial slice cuts
keep its three matter spins together, so the full-state spatial entropy is
unchanged. Physical local matter observables require the inverse basis map;
mode occupancies must not be relabeled individual physical matter Sz values.
Gauge Sz and total matter number are unchanged. Generic finite observable
and winding-measurement drivers map these states back to the physical basis.
Checkpoint metadata explicitly names the basis to prevent mixing conventions.

`scripts/winding_qn_dmrg.jl` is the direct-constraint optimization driver.
Every update remains in the number/winding block. It records measured diagonal
loop means/variances, which represent the original matter-dressed microscopic
loop by the verified basis identity, and validates energy variance, sweep,
truncation and entropy convergence. The exact projector-penalty driver remains
an independently validated fallback and cross-check. No direct-QN optimized
state is claimed until the library tests and sector ED comparisons pass.

For an infinite cylinder, periodically repeating onsite weights is not yet
proven to fix a single winding eigenvalue: it can instead fix a product of
several translated generators. No such periodic charge is silently treated
as a definite global flux/MES. Resolve this before extending the direct-QN
construction to infinite sector calculations.

A proposed checkpoint bridge avoids restarting wider states or forming their
full wavefunctions. `src/winding_checkpoint_bridge.jl` uses a library-decomposed
local B† MPO and a three-state virtual charge accumulator to implement the exact
operator P_q B† between original U1 site indices and the new U1 x Z3 site
indices. The library performs MPO application and compression. Projection
weight, number preservation, cutoff and source sector must be recorded.
All 36 checks in `test/winding_checkpoint_bridge.jl` passed, including independent physical-basis projector weights, state overlaps, number preservation and energy equivalence for both families and all three winding charges. Original states remain immutable; conversion produces a new branch.

The first direct optimizer attempt passed all 42 physical-library basis tests,
then reached E=-2.6131259237 in a restricted N2,q0 variational branch and failed
with KrylovKit's LAPACK STEGR `LAPACKException(22)` during a local solve. This
is not a failed Hamiltonian/basis validation or a ground-energy result. Its
latest valid checkpoint is preserved by an immutable byte-copy snapshot.
The next continuation uses a smaller Krylov subspace (12), more restarts (30),
and a recorded density-matrix noise schedule [1e-5,1e-6,1e-7,0] to expand the
allowed QN variational support. The requested local tolerance 1e-13 and final
convergence criteria remain unchanged. All updates still preserve physical
number and the target winding QN. The known exact-sector physical states also
provide a checkpoint-bridge initialization route if the product/random QN
initializer remains trapped.

The resumed direct N2,q0 run now reaches the independent ED minimum
-2.73463474709466, full-state S=0.50662555964812, energy variance within
roundoff and winding-purity error zero. No requested tolerance was loosened.
The checkpoint bridge passed all 36 tests after correcting signed number-flow
labels for library MPO links pointing In. The correction concerns the bridge
operator, not the microscopic Hamiltonian or earlier energy optimizations.

Known small-sector ED energies are now an explicit optimizer acceptance gate,
so a stationary excited branch cannot be labeled a passed ground candidate.
With the recorded startup noise schedule, truncation acceptance uses the last
four zero-noise sweeps at the same 1e-10 threshold; all noisy and zero-noise
truncation errors remain saved. The physical-filling optimization has not yet
been certified. Additional-cycle ED resolves different subblocks, but no
saved-state extra-cycle charge is inferred from its energy alone.

The bridged armchair N2,q0 continuation reproduced its exact ED ground energy and zero energy/loop variance at every stage, but repeated density-matrix noise rotated its four-dimensional ground space. Its entropy changed by 0.0010–0.0014 between stages, so the unchanged 1e-5 entropy criterion correctly rejected convergence. This is preserved in `direct_winding_qn_bridged_fixtures_and_physical_sectors.log` and the individual result JSON, rather than labeled a converged entropy. Continuation now uses zero noise on that saved state. For subsequent new starts, expansion noise is applied only in the first stage; later stages have zero noise. This separates initialization from degenerate-state stability without relaxing acceptance thresholds. A stable chosen state is still not a proved MES.

An independent eight-state matrix identity also gives a real closed form for the rotated Hamiltonian, recorded in `real_charge_hamiltonian_audit.json`. The proposed library MPO implementation is not used in production until its equivalence checks pass. All running optimizations retain their loaded implementation.
