# Exact microscopic winding-sector optimization

Priority: validate a sector-preserving variational construction before new
MES, topology, or deconfinement studies. The operators here are exact CGS
winding symmetries. Identifying their eigenvalues with emergent topological
flux is a separate physical question; local cycle purity is not that proof.

For a recorded winding-one exponent vector p, every vertex obeys
p_1+p_2+p_3=0 mod 3. Write p_i=c+q(i-1) mod 3 at each vertex. The operator
acts on all three local matter spins by their cyclic permutation q and the
common physical number rotation c, and on each shared gauge spin by its
recorded number rotation p_e. Each gauge edge is acted on once. The local
W intertwining identity establishes commutation with the exact exchange
Hamiltonian. U^3=I and U preserves N_up. The explicit endpoint actions are
recorded in the geometry cycle JSON files.

The physical individual-spin MPS QN implements U(1), while a winding U contains
three-spin matter permutations and is not a product of individual-spin
on-site charge labels in this basis. A new local charge basis is not assumed.
Use instead

P_q=(I+omega^(-q)U+omega^(-2q)U^2)/3,
H_q=H+lambda(I-P_q).

This is a Hermitian orthogonal projector penalty. H_q restricted to sector q
is exactly H, with no rescaling or added energy. If B=sum of exchange-term
norms=3 sqrt(3) Nv, then ||H||<=B. Choosing lambda=2B+1 places every undesired
sector strictly above every state in the target sector. This conservative
choice has conditioning costs, which must be checked rather than silently
weakening it. The sum-of-MPO DMRG interface avoids a truncated MPO addition.
No custom tensor-network algorithm or giant circumference supersite is used.

Individual DMRG updates can still contain sector mixtures. This construction
isolates the exact target eigenspace in the optimization problem; finite-bond
stationarity, leakage, and energy variance must be validated afterwards.
It is not an assertion that every truncated intermediate tensor is pure.
If purity fails, refine from that checkpoint or implement stronger constrained
updates; do not describe an impure final state as sector-preserving success.

Independent ED in complete N_up=2 blocks on zigzag L2,w1 (19 spins) and
armchair L2,w1 (20 spins) validates all three projectors, Hamiltonian
commutation, U^3, Hermiticity, and unchanged within-sector Hamiltonian.
These inexpensive blocks test the algebra and penalty, not ground filling.
`results/sector_penalty_ed_audit.json` preserves errors, target energies,
sector dimensions, fixed filling and provenance. Julia MPO validation passed all 36 new checks at the stated tolerances;
24 translated circumference checks and three local-loop checks also passed.
The separately optimized N_up=2 q0/q1 fixtures now pass for both geometries;
physical ground-number sector optimization is the next validation.

`scripts/sector_preserving_dmrg.jl` records N_up, physical S^z, filling, exact
cycle definition, penalty, sweeps, errors, entropy/Schmidt spectrum, loop mean
and variance, and physical-H energy variance. It checkpoints every two sweeps,
keeps completion files immutable, and requires 1e-10 loop purity, 1e-8 physical
energy variance, 1e-9 final energy drift, 1e-5 entropy drift and 1e-10 truncation
before reporting convergence in the fixed number/loop sector. These criteria
do not certify global filling selection, MES, length convergence or 2D order.

Commands (after MPO tests pass):

```sh
OPENBLAS_NUM_THREADS=1 python scripts/sector_penalty_ed_audit.py
julia --project=. test/cgs_sector_penalty.jl
julia --project=. scripts/sector_preserving_dmrg.jl zigzag 2 1 2 0 64
julia --project=. scripts/sector_preserving_dmrg.jl zigzag 2 1 2 1 64
julia --project=. scripts/sector_preserving_dmrg.jl armchair 2 1 2 0 64
julia --project=. scripts/sector_preserving_dmrg.jl armchair 2 1 2 1 64
```

The low-number MPS fixture energies must agree with independent ED. Once
validated, optimize the physical candidate number sectors (zigzag N9/10,
armchair N8/12 and nearby sectors as needed) separately. Preserve all previous
unconstrained optimizations. Finite MPS do not provide a unique infinite
transfer length by joining their end bonds arbitrarily.

Infinite extension is not yet implemented. Penalizing every translated winding
loop can impose additional contractible-sector constraints, so it must not be
silently identified with fixing a single global flux. That extension requires
an explicit symmetry/sector specification and its own validation.

A new no-diagonalization audit verifies antiunitary covariance in paired number
blocks. Particle-hole maps (N,q) to (Ns-N,sum_e p_e-q mod 3). Both recorded
winding loops have sum_e p_e=0 mod 3, so q maps to -q. This equates paired
sector energies, and does not prove half filling minimizes energy. Errors and
fixed-sector conventions are preserved in `results/sector_particle_hole_audit.json`.

An exact local matter basis change now supplies a promising direct U1 x Z3
implementation. All eight matter states and individual two-dimensional sites
are retained. Local algebra and independent full-block ED passed; library
validation is pending. See `results/direct_winding_qn_strategy.md`. This is
the preferred continuation if its production checks pass; the projector
penalty remains an exact fallback and independent cross-check.
