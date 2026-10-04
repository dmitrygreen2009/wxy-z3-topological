# Exact microscopic WXY Z3 CGS numerics

The model uses three local matter spin-1/2 degrees of freedom per honeycomb
vertex and one shared gauge spin per physical edge. All calculations retain
matter and gauge spins. Production MPS calculations use ITensors.jl /
ITensorMPS.jl; optional VUMPS uses ITensorInfiniteMPS.jl.

**Torus spectrum audit:** the supplied eight-level list omits degenerate
partners. Independent Julia block Lanczos, Python ED, and separate exact
parity blocks establish four ground states for the stated full N_up=9 block.
See [the complete hard-gate audit](results/torus_hard_gate_audit.md).
The corrected first ten levels, including four ground states and four states
at -7.296046028855, are now the regression benchmark. Cylinder scaling resumed
after the hard gate was resolved; the audit remains the historical record.

## Reproduce

Use Julia 1.10.10 and the committed Project/Manifest. Python ED is independent.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python scripts/ed.py
julia --project=. test/runtests.jl
julia --project=. scripts/validate.jl
julia --project=. scripts/batch.jl
julia --project=. scripts/refined_batch.jl
julia --project=. scripts/infinite.jl zigzag 1 32
```

`scripts/cylinders.jl FAMILY LENGTH WIDTH CHI` runs an individual finite case.
Append `star` to use the improved ordering within each spatial slice. The
spatial cut contains exactly the same physical sites under either ordering;
each site remains an individual spin-1/2. Star ordering interleaves nearby
matter and gauge spins to reduce entanglement at intermediate MPS bonds.
`scripts/continue.jl CHECKPOINT.jls CHI` increases its bond dimension.
`scripts/observables.jl CHECKPOINT.jls ...` measures full-state entropies and
connected matter correlations. Checkpoints are local, git-ignored Julia
serialization files, requiring the pinned environment. JSON results and logs
are preserved in `results/`.

## Lattice conventions

Nearest-neighbor honeycomb distance is one. Let triangular Bravais vectors
`a1=(sqrt(3)/2,3/2)`, `a2=(-sqrt(3)/2,3/2)`, and displacement `B-A=(0,1)`.
A(x,y) has legs to B(x,y), B(x-1,y), B(x,y-1). The corresponding gauge edge
is keyed by the A coordinate and its leg; its B endpoint uses that same leg.

Zigzag circumference is `w*a2`, length `sqrt(3)*w`; axial slice is `t=x`.
Armchair circumference is `w*(a1+a2)`, length `3*w`; axial slice is `t=x-y`.
Open ends retain dangling gauge spins so every star has the exact three legs.
These add boundary spins to the bulk count `4.5*Nv`. Periodic parallel edges
are retained, including at width one. Finite entropies are natural logarithms
of the full-state Schmidt probabilities at the central transverse cut.

The MPS contains individual spin sites, sorted by axial position. No whole
circumference is made into a supersite. Matter spins of a vertex remain on one
side of the spatial cut; boundary-crossing gauge spins are assigned by their
edge midpoint, with a fixed small axial tie break. See the explicit
[four-star incidence table](results/torus_incidence.md).

Finite-MPS transfer spectra are not defined by arbitrary identification of
unrelated bond spaces. Infinite calculations use the library's TransferMatrix;
their correlation length is reported in two-slice cells and in axial slices.
The corresponding axial distance per slice is 3/2 (zigzag) or sqrt(3)/2
(armchair).

## Environment and regression tests

The recorded environment is Julia 1.10.10 and Python 3.14.0; complete package
versions are in `results/environment.json`, Julia's Project/Manifest, and
`requirements.txt`. Use `JULIA_DEPOT_PATH` and a Python virtual environment
outside tracked scientific output. A new machine can run:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/ed.py
julia --project=. test/runtests.jl
julia --project=. test/entanglement.jl
julia --project=. test/transfer.jl
julia --project=. test/checkpoints.jl
julia --project=. scripts/torus_independent.jl
julia --project=. scripts/validate.jl
```

Normal GitHub CI runs inexpensive ED, geometry/Hermiticity, entropy and
analytic transfer tests. The full independent block-solver audit and ITensor
benchmark DMRG run on manual workflow dispatch. Cylinder DMRG/VUMPS is never
part of normal CI. The corrected first ten torus levels have multiplicities
4,1,1,4, including four ground states; the audit remains unchanged.

## Complete geometry files

`geometry/` contains machine-readable benchmark and cylinder graphs, with
vertices, local matter species, physical edge-spin IDs, both endpoint leg
indices, A/B conventions, periodic translations, and spatial partitions.
Infinite files store translated interaction pairs and the two-slice unit
cell. Regenerate them with:

```sh
julia --project=. scripts/export_geometry.jl
```

New finite production cases export their own geometry before optimization.
The microscopic Hamiltonian is
`H = -sum(v,a,i)[W[a,i]*Sminus(matter[v,a])*Splus(edge[v,i]) + conjugate(W[a,i])*Splus(matter[v,a])*Sminus(edge[v,i])]`,
with `W[a,i]=omega^((a-1)*(i-1))/sqrt(3)` and `omega=exp(2*pi*im/3)`.
There are three local matter spins per vertex and one shared physical gauge
spin per edge. No effective clock or gauge-only entropy model is used.

## Checkpoint and resume

Production version 3 stores local checkpoints under `results/checkpoints/`,
with geometry, circumference/width, physical size, bond dimension, sector,
seed, and run version in names. Small JSON manifests contain SHA-256 hashes
and point to the payloads. Periodic writes are atomic and retain a previous
payload/manifest pair. Checkpointing occurs every two DMRG sweeps and every
two VUMPS iterations. Completed result JSON files point to their descriptive
checkpoint. Legacy checkpoint names remain readable through the resolver.

```sh
julia --project=. scripts/cylinders.jl zigzag 6 2 256 star
julia --project=. scripts/cylinders.jl armchair 6 2 256 star
julia --project=. scripts/continue.jl results/zigzag_L6_w2_chi256_star.jls 512
julia --project=. scripts/resume_checkpoint.jl --latest zigzag 6 2 512
julia --project=. scripts/resume_checkpoint.jl results/checkpoints/NAME.jls
julia --project=. scripts/infinite_qn.jl zigzag 1 64
julia --project=. scripts/infinite_from_finite.jl results/zigzag_L6_w2_chi256_star.jls
julia --project=. scripts/continue_infinite_audited.jl results/infinite_zigzag_w1_chi32_qn_star_active.jls
```

Do not restart an interrupted large job before checking for the newest valid
checkpoint. The resume command verifies checksum, size and normalization.
Finite correlation length is left null rather than identifying unrelated
finite bond spaces. Infinite measurements save complex leading transfer
values, residuals, canonical/isometry checks, and xi in two-slice cells and
axial slices. Values below numerical transfer resolution are flagged.

## Analysis, fits and provenance

```sh
OPENBLAS_NUM_THREADS=1 python scripts/scientific_summary.py
OPENBLAS_NUM_THREADS=1 python scripts/fit_entropy.py
MPLCONFIGDIR=/tmp/wxy-mpl python scripts/analyze.py
```

Raw JSON/CSV stays separate from `entropy_fits.json`. Each fit records its
circumferences, slope, gamma, OLS scatter uncertainty, residuals, comparison
to ln(3), and the fit with its smallest circumference removed. A two-point
fit has no residual degrees of freedom and therefore no scatter uncertainty.
These uncertainties exclude systematic length, chi, circumference and flux
sector effects. Unstable intercepts do not establish or rule out topology.

Production audit records include seed, solver/settings, initialization,
conserved quantities, source fingerprints, git revision, dirty-code status,
runtime and numerical errors. Older runs predate complete launch metadata;
these gaps must be identified rather than inventing missing provenance.
Resource estimates are recorded before increases. The conservative per-worker
memory budget is half of system RAM; existing smaller checkpoints survive a
rejected larger request. Large payloads and caches remain local. Validated
milestones are committed and pushed to the current branch at natural
checkpoints. A push failure is recorded and never discards local work.

Historical iteration data can be recovered with
`python scripts/archive_legacy_logs.py`. Its CSV preserves log line identities
and runtime lower bounds; the accompanying audit states missing provenance.
New jobs capture source fingerprints once at process launch.

After `scripts/extend_study.jl` finishes, `julia --project=.
scripts/additional_convergence.jl` checks a common bond dimension of 512 for
widths one through four and longer wide cylinders. Resource estimates can skip
a point while preserving all completed smaller states; a resource skip is
recorded separately from numerical convergence. Independent seed checks use
`julia --project=. scripts/multiple_seeds.jl`.
