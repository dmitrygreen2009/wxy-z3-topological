# Microscopic circumference-operator measurement validation

The new `src/cgs_operator_mpo.jl` constructs exactly the matter permutations and
physical shared-gauge rotations in the machine-readable winding-one cycles.
It uses ITensorMPS local three-spin MPO decomposition and unit-bond identities,
with no dense whole-circumference operator and no custom contraction algorithm.
It requires the audited individual-spin star ordering.

`test/cgs_operator_mpo.jl` passed all 24 checks for zigzag and armchair widths
one and two: normalized action and complex overlap against the independent
`cgs_cycle_gates` implementation on random finite states, operator bond dimension
at most four, and exact product-state expectations at offsets zero, one slice,
and one complete infinite unit cell. Tolerances remain 1e-11 or tighter.

Two implementation defects were resolved before production use. CelledVector
range indexing rejected translated support even though scalar indexing supports
it. The scalar-index construction fixes this. The library's canonicalized
identity-product MPO carries compensating scalar phases; replacing individual
canonicalized nodes loses those phases. Building local identity tensors on the
library's unit-bond scaffold preserves the operator phase. The failed test logs
are retained alongside the passing validation log; no stored production state
or previously validated numerical result was changed.

The charge-three matter correlation operator is the exact product of three S+
operators on one star and three S- operators on another. Independent coherent
product expectations are 1/64; an all-down product gives zero. A Hadamard acting
on |Dn> yields the minus coherent state, whose triplet one-point expectation is
-1/8, while its two-triplet expectation remains 1/64. These signs are fixed
analytically rather than by relaxing any tolerance. All five new correlation
checks passed at tolerance 1e-12, including support beyond one infinite cell.

`scripts/infinite_circumference_cycles.jl` measures both slice-origin winding
operators and charge-three correlations at separations one, two, and four
infinite cells. It saves complex expectations, charge probabilities, source
payload checksums, inherited canonical/solver errors, and provenance. These are
observables of the saved variational state. Neither a pure loop charge nor a
short variational correlation length alone certifies a minimally entangled
state, a spectral gap, or the two-dimensional phase.

## First saved-state measurement (provisional)

The zigzag width-two finite-seed state at cap 32 (actual bond dimension 16) was
recanonicalized from the same AL-defined state because its stored center
inconsistency was 2.83e-4. Canonical error is now 6.36e-12; its inherited
variational solver residual remains 2.83e-4. No optimization was repeated.
Its two winding operators have charge weights approximately
(0.356, 0.322, 0.322) and (0.363, 0.319, 0.319), rather than a pure charge.

The physical, CGS-invariant matter-triplet expectations have magnitudes about
0.0432 and 0.0460. At separation four cells the raw correlations are about
0.00187 and 0.00212, while connected parts are near floating-point roundoff.
Thus this particular unrestricted variational ansatz breaks physical U(1),
with a disconnected ordered plateau and short connected correlations.
This is not a claim of spontaneous order in the converged cylinder or 2D
model. Its finite bond dimension, residual, initializer branch, and sector
weights still require convergence checks. In particular its short transfer
length is not evidence of a gapped topological phase. The ongoing wider-bond
run is retained for comparison, rather than restarted.

Raw expectations remain in the `_circumference_cycles.json` result.
`scripts/summarize_infinite_microscopic_observables.py` generates a separate
provisional summary without promoting these observations to phase conclusions.

## Contractible charges and number-sector follow-up

The lifted hexagon exporter records 450 closed walks with zero unwrapped
winding, exact physical edge IDs, and local endpoint transformations. Twelve
new tests passed against independently applied physical-spin gates for both
wrappings and widths one and two, including infinite translated support.
The width-two saved state's contractible hexagon charge weights are also
mixed, approximately (0.342, 0.329, 0.329) and (0.335, 0.333, 0.333).
These measurements do not assume that a microscopic ground sector must have
charge zero.

At armchair L4 width one, the refined charged Nup18 candidate has energy
-15.345428885380793 and full-state entropy 2.183249996832878, compared with
-15.34530743538934 and 2.1949203117244993 for the Nup19 candidate. Their last
sweep energy drifts are about 2e-13, but truncation errors 1.30e-8 and 1.97e-8
still exceed the strict 1e-8 acceptance target. No tolerance was loosened.
The charged state's first winding loop has q1 while the neutral state's has
q0; other measured winding loops have q0. Thus both number and flux content
change between these variational candidates. The number difference in the
central physical interval [1,3) is -0.578 of the total -1; this L4 measurement
cannot establish a boundary-only or bulk density effect. Independent charged
initialization of an unrestricted infinite width-two armchair state now tests
the filling assumption directly. Exact neutral and charged checkpoints are
preserved in local byte copies with committed small snapshot manifests.

Geometry records now distinguish exact physical axial coordinates from the
legacy 0.01/0.1 ordering keys. All physical cuts are unchanged, and all 41
finite and 16 infinite static geometry checks pass. No validated energy,
Schmidt spectrum, or production optimization was repeated for this metadata
clarification.

## Separated loop products and transfer-length interpretation

Twelve further checks validate disjoint products U(first) U(second)^dagger
against independent finite physical-spin gates, normalization, and translated
infinite product expectations for both wrappings at widths one and two.
The adjoint uses the exact Abelian representation (c,q,p) -> (-c,-q,-p) mod 3.
No whole-circumference dense tensor or custom contraction kernel is introduced.

New charge-three connected contractions on the saved, consistently canonical
zigzag width-one QN cap128 state are only about 1e-10 at one/two cells, 1e-11
at four cells, and 2e-15 at eight cells. These are observables of the approximate
state; their small amplitudes are not a physical-ground-state correlation error
bound, and the inherited solver residual remains 3.26e-5. They do not demonstrate
long charge-three correlations even though the full transfer length is about
17 axial slices. Additional winding-loop pair contractions test whether that
full transfer mode instead couples to narrow-quotient sector fluctuations.
Existing validated one-point contractions are reused rather than repeated.
