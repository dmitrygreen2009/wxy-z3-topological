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
