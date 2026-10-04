# 18-spin torus hard-gate audit

## Resolution

The supplied first-eight list is **not the complete low spectrum** of the
stated microscopic Hamiltonian in the full 48,620-dimensional N_up=9 sector.
It undercounts the ground multiplicity and two excited multiplicities.
Independent Julia block Lanczos, separately constructed exact parity blocks,
and Python ED agree on four ground states. This meets hard-gate condition (2),
an independent numerical demonstration that the supplied spectrum is wrong
for the stated full sector. It does not reproduce the supplied list as a
valid complete spectrum, and no tolerance was loosened to declare it passed.

| Index | Full-sector result | Supplied list |
|---|---:|---:|
|0|-7.364534839426|-7.364534839426|
|1|-7.364534839426|-7.364534839426|
|2|-7.364534839426|-7.342883075944|
|3|-7.364534839426|-7.337288472473|
|4|-7.342883075944|-7.296046028855|
|5|-7.337288472473|-7.296046028855|
|6|-7.296046028855|-7.289970027036|
|7|-7.296046028855|-7.289970027036|

The first fourteen levels have multiplicities **4, 1, 1, 4, 4** at
`-7.364534839426`, `-7.342883075944`, `-7.337288472473`,
`-7.296046028855`, `-7.289970027036`, respectively.
The maximum discrepancy between the *ordered full first eight* and the
supplied list is 0.046837047089, despite agreement on E0.

## Independent eigenstate check

`scripts/torus_independent.jl` does not import the Python Hamiltonian, Python
eigenvectors, or ITensor model constructor. It constructs the sparse operator
directly from an endpoint-and-leg edge table; enumerates all 2^18 kets by
popcount, with reversed bit order; and uses eight-vector **block Lanczos**
rather than scalar ARPACK. It converges sixteen eigenpairs.

Its four lowest energies are:

```
-7.364534839425768
-7.364534839425765
-7.364534839425747
-7.3645348394257395
```

Their residual norms are below 8.0e-14. The sixteen-vector Gram-matrix error
`||V†V-I||` is 1.21e-13, and all sixteen residuals are below 7.9e-13.
Thus these are four independent states, not copies of one or two vectors.
The first excited energy is separated by about 0.02165, far beyond these
numerical errors. An orthonormal four-dimensional subspace at E0 also gives
a min-max certificate that a two-ground-state spectrum with its third level
0.02165 higher cannot describe this operator.

Raw results: `torus_independent.json`, `ed.json`.

## Geometry and shared-edge audit

Ordinary periodic honeycomb translations preserve A/B sublattices. Four
vertices mean two primitive cells. Every index-two translation subgroup of
Z² is the kernel of one of the three parity maps
`x mod 2`, `y mod 2`, `(x+y) mod 2`. All three were constructed explicitly
and diagonalized. Their full spectra, including multiplicities, agree.

The local-leg tables, in vertex order A0, B0, A1, B1, are:

| Parity map | A0 | B0 | A1 | B1 |
|---|---|---|---|---|
|x|[0,1,2]|[0,4,2]|[3,4,5]|[3,1,5]|
|y|[0,1,2]|[0,1,5]|[3,4,5]|[3,4,2]|
|x+y|[0,1,2]|[0,4,5]|[3,4,5]|[3,1,2]|

Indices here are zero-based physical gauge-edge IDs. Every edge appears at
exactly two vertices, and each star has three distinct incident edge IDs.
There are twelve matter spins and **six gauge spins**, not twelve independent
gauge endpoint copies. The gauge spin assigned to either endpoint of edge e
is the same physical site 12+e. Hermiticity is checked explicitly.

Parallel edges are retained. A connected cubic bipartite multigraph with
two A and two B vertices has multiplicity matrix
`[[k,3-k],[3-k,k]]`; connectedness forces k=1 or 2. These are the same graph
up to relabeling. k=0 or 3 would be disconnected and cannot be a quotient of
the connected honeycomb. K4 is non-bipartite and cannot be this cell.

See `torus_incidence.md` and the endpoint tables in
`torus_convention_audit.json`.

## Exhaustive endpoint orientation audit

For each star independently, all six leg permutations and both W/conj(W)
choices were considered: 12^4 patterns for each of three translation cells,
or **62,208 labeled cases**. The complete first eight for every pattern are
recorded in `torus_all_endpoint_conventions.csv.gz`. None matches the supplied
complete list.

Equivalent cases were reduced using an explicit local-unitary certificate,
not assumed equivalent or individually rerun 62,208 times. For zero-based
indices, `W[a,i]=omega^(a*i)/sqrt(3)`. Every permutation of three columns is
affine, `p(i)=s*i+t mod 3`, with s=+1 or -1. Therefore

```
W[a,p(i)] = omega^(a*t) W[s*a mod 3,i].
conj(W[a,i]) = W[-a mod 3,i].
```

These are permutations and phase rotations of the **local matter spins**.
They are implemented by matter-label permutation unitaries and single-spin
z rotations. They preserve N_up and leave every shared gauge spin untouched,
so they preserve the entire spectrum and all multiplicities even when each
endpoint chooses a different convention. The twelve monomial identities were
checked numerically below 1e-14. Arbitrary matter-label permutations are also
basis relabelings. W is symmetric, so transposing W adds no new case.

No extra projection, flux restriction, boundary twist, or gauge approximation
was introduced. An independently specified twisted Hamiltonian or extra
sector restriction would be a different benchmark.

## Symmetry and eigensolver audit

Three independent exact Z3 cycle operators were verified on the microscopic
torus, including their required matter rotations and permutations:
`max|[H,U]| < 3.6e-16`, `max|U³-I| < 1.1e-15`.
They arise from the actual periodic incidence, and occur in all equivalent
geometries; changing an endpoint leg convention cannot remove them.
Their ground-space leakage is below 2.1e-13. These checks do not infer topology.

There is also an exact Z2 symmetry R: flip all eighteen spins, then swap
matter labels 2 and 3 at every star. R preserves N_up=9 and commutes with H.
Its + and - parity blocks are disjoint 24,310-dimensional orthogonal sectors.
**Each block has two ground states at E0**, establishing four in their union.
The -7.342883075944 singlet lies in the + block; the -7.337288472473 singlet
lies in the - block. Hence the supplied list is not simply one parity sector.

Scalar ARPACK was checked on a declared grid of tolerances, Krylov dimensions,
and seeds. It returned one, two, or four ground vectors depending on settings.
With `k=8, ncv=20, tol=1e-10`, and the documented seed-one complex start,
it reproduces **all eight supplied numbers** to 1e-9 while returning only two
of the four ground partners. Its returned eigenvectors have small residuals,
but this does not certify completeness of a degenerate eigenspace.
This demonstrates a mechanism for exactly the discrepancy. It does not prove
which solver originally produced the supplied benchmark.

Raw parity and scalar-solver results are in `torus_solver_audit.json`.
The corrected regression checks all fourteen low levels, residuals, and
orthogonality. CI additionally runs the independent Julia and convention audits.

## Current scope

Cylinder, VUMPS, and DMRG scaling remain stopped as requested. Existing larger
results are exploratory and establish no scientific conclusion about Z3
topological order. The audit resolves the torus discrepancy for the exact
Hamiltonian and full fixed-number sector specified in the request.
