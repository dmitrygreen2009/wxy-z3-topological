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
  of the full-state spatial Schmidt entropy. Execution is pending.

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
