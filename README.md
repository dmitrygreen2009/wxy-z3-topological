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

**Current scientific status:** topology, a 2D gap, and physical U(1) order remain unresolved. Read [the strategy re-audit](results/scientific_strategy_reaudit.md) before interpreting any apparent intercept or transfer length. Static geometry/cut checks passed; a projected initializer was also found to change its CGS sector during DMRG. Existing numerical states are preserved, and none of the current gamma fits is accepted as a thermodynamic topological estimate.

## Reproduce

Use Julia 1.10.10 and the committed Project/Manifest. Python ED is independent.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
python -m pip install -r requirements.lock.txt
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
.venv/bin/python -m pip install -r requirements.lock.txt
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

To audit seed-dependent loop trapping, use `julia --project=.
scripts/cylinder_loops.jl RESULT.jls` and `scripts/check_finite.jl RESULT.jls`.
`test/cgs_projection.jl` validates exact microscopic cycle projection.
`julia --project=. scripts/projected_zigzag.jl 4 6 10` tests trivial-cycle
initializations for zigzag width one through chi=512. This selected sector
must be compared against unrestricted runs; no global optimality at other
widths or lengths is assumed.

Infinite entropy consistency can be audited independently with
`julia --project=. test/infinite_entropy.jl` and
`julia --project=. scripts/recanonicalize_infinite.jl RESULT.jls`.
Recanonicalization uses the official transfer-fixed-point routines, records
inherited variational errors, and skips nonunique dominant fixed points to
preserve sector weights. It does not certify ground-state convergence.
`python scripts/benchmark_number_sectors.py` scans all particle-number sectors
of the small exact benchmarks, independently of the torus multiplicity test.

Additional independent finite checks:
```sh
python scripts/benchmark_number_sectors.py
python scripts/flatband_audit.py
python scripts/small_cylinders_ed.py
julia --project=. scripts/cylinder_block_ed.jl zigzag armchair
```
The last command uses eight-vector block Lanczos on the explicit L2 width-one
geometry records to retain low-energy multiplicities. Scalar-ARPACK results
certify their returned eigenpairs, rather than complete degeneracy counts.
Small-cylinder ED eigenvector entropy can depend on the chosen degenerate
ground-state combination. These open-cylinder spectra do not establish a
thermodynamic phase.

For a finite-to-infinite warm start, the library canonicalization can reject
a positive imaginary transfer eigenvalue component above 1e-15 even when
the eigensolver converges at 1e-14. The driver retries seeded library Krylov
starts on the same cached fitted state, preserving all tolerances. Every
attempt and exact error is saved in the fit's canonicalization JSON.

Refresh the additional full-state spectral, infinite-convergence, and plotting
summaries from saved small results:
```sh
python scripts/entanglement_spectra.py
python scripts/analyze_infinite.py
MPLCONFIGDIR=/tmp/wxy-mpl python scripts/plot_study.py
```
`infinite_analysis.json` separates raw measurements, selected canonical
remeasurements, bond-dimension changes, and circumference fits. The two-point
finite-entanglement effective-central-charge slope is descriptive and cannot
establish a critical bulk phase. Missing circumference fits indicate too few
widths; no synthetic data are inserted.

Number-sector selection is part of the ground-state audit. The complete tiny
open-cylinder scan finds an armchair global minimum outside half filling;
see `results/number_sector_cylinder_audit.md`. Existing half-filled armchair
points describe fixed-number-sector states. To run a specified sector:
```sh
julia --project=. -e 'include("scripts/cylinders.jl"); run_cylinder("armchair",2,1,128;ordering="star",nup=8)'
python scripts/benchmark_number_sectors.py --cylinders
```
Explicit-sector result names include `_Nup8` and checkpoints verify the MPS
quantum-number flux. Continuation and explicit checkpoint resume retain the
sector. The default `--latest` resume chooses the closest half-filled sector;
use an explicit checkpoint path for a different sector.

Independent exact symmetry audits:
```sh
python scripts/parallel_edge_symmetry.py
julia --project=. test/particle_hole.jl
julia --project=. scripts/particle_hole_audit.jl results/armchair_L4_w1_chi256_star_charge_sectors_global_audit.json
```
`particle_hole_mps` implements the exact antiunitary product of physical spin
flips and complex conjugation using library tensors, preserves QNs, and maps
N_up to N-N_up without changing energy or spatial entropy. This supplies
partner-sector starting states and exposes independent solver discrepancies.
The narrow zigzag parallel-edge graph has additional local S3 symmetry; see
`results/infinite_sector_audit.md`. Its mixed-cycle infinite measurements
must not be pooled with the finite projected q=0 states as one sector.

The star ordering sorts the declared physical left/right spatial sets before
its local star order, so odd as well as even finite lengths have the correct
individual-site entropy partition. Existing even-length production orders
are unchanged. The finite-to-infinite bulk-slice fitter requires an even
length of at least four; this is checked explicitly. A permanent regression
compares both orderings and physical spatial membership on 24 geometries.
Checksummed manifests are verified on state loads and continuations; explicit
`.previous` payloads retain their own state instead of following a newer
manifest. Legacy unmanifested fitting caches remain explicitly identifiable.

Completed finite checkpoints now have unique completion identifiers; rolling latest checkpoints still rotate a previous payload. `scripts/star_trial_bound.py` independently validates the isolated-A-star variational trial state. Its conservative two-slice bound is `E0 <= -4.8061842 * width` for either wrapping. A normalized, canonically consistent variational state above that trial energy cannot be the global ground state, even if its projected residual is tiny. Inconsistent stored center tensors do not justify that exclusion from the raw expectation; fixed-density exclusions additionally require a density-compatible trial. `test/variational_bound.jl` additionally checks the disjoint-dimer trial state directly with ITensor.

Re-audit geometry and cycle winding without rerunning the numerical solvers:

```sh
python scripts/strategy_geometry_audit.py
OPENBLAS_NUM_THREADS=1 python scripts/cgs_general_audit.py
julia --project=. scripts/noncontractible_cycles.jl results/zigzag_L4_w2_chi256_star.jls
```

The last command measures new circumference-cycle charges and density profiles on an existing state. CGS projection labels specify initialization; verify actual charges after optimization. A winding-operator eigenstate has not automatically been identified as a topological MES. Future finite-to-infinite fit caches carry source and payload fingerprints; legacy caches with unknown input fingerprints remain explicitly labeled reusable trial states with incomplete initializer provenance.

### Audited circumference charges

For a saved infinite state in the explicit `star` ordering, measure winding-one
microscopic CGS loops with:

```sh
julia --project=. test/cgs_operator_mpo.jl
julia --project=. scripts/infinite_circumference_cycles.jl RESULT.jls
```

The operator is factored into local three-spin library MPOs, so its construction
does not allocate a dense whole-circumference tensor. Results include complex
loop expectations and three charge probabilities,
source checksums, and inherited convergence errors. Use a consistent canonical
state for these contractions; `scripts/recanonicalize_infinite.jl` rebuilds the
same AL-defined state without optimizing it when stored centers are inconsistent.
See `results/circumference_operator_validation.md` for the independent checks.
Loop charge alone does not identify a minimally entangled state or certify a
2D phase, and unconverged variational correlations remain provisional.

When another running driver may reuse a completed result alias, preserve its
exact local checkpoint and small result manifest first:

```sh
python scripts/preserve_result_snapshot.py results/RESULT.json descriptive_label
```

This verifies the input checksum and copies bytes without repeating any
calculation. Large payloads remain git-ignored; snapshot JSON records are kept
under `results/snapshots/`.

Contractible hexagons are stored separately, with explicitly closed unwrapped
walks and zero winding. Their operator and saved-state measurements use:

```sh
python scripts/export_contractible_cycles.py
julia --project=. test/contractible_operator_mpo.jl
julia --project=. scripts/infinite_contractible_cycles.jl RESULT.jls
python scripts/summarize_infinite_microscopic_observables.py
```

Both contractible and circumference charge probabilities need convergence
checks. The two kinds of cycle have separate machine-readable definitions
under `geometry/cgs_plaquettes/` and `geometry/cgs_cycles/`. Physical axial
coordinates are recorded explicitly; the ordering tie-break keys are not
physical distances. `scripts/compare_number_density_profiles.py` compares
saved neutral/charged density measurements using these physical coordinates.

Legacy in-flight drivers may still emit the older ordering-key description.
At their next completed checkpoint, annotate those metadata fields and rerun
only the static geometry checks before exporting cycle definitions:

```sh
python scripts/normalize_geometry_metadata.py
python scripts/strategy_geometry_audit.py
python scripts/export_contractible_cycles.py
```

This changes descriptive metadata only; physical edges, endpoint legs, spin
ordering, and spatial cuts are checked independently and remain preserved.

To compare the full transfer length with new physical or narrow-sector
correlations while reusing previously measured one-point data:

```sh
julia --project=. test/cgs_pair_operator_mpo.jl
julia --project=. scripts/infinite_physical_correlations.jl RESULT.jls RESULT_observables.json
julia --project=. scripts/infinite_sector_correlations.jl RESULT.jls RESULT_observables.json
```

The last command is specifically for the zigzag width-one quotient. Its
winding operators have unusually local support, so these sector correlations
must not be interpreted as evidence for the two-dimensional phase.

An independent infinite-state initializer tests fixed filling without relying
on the product-state expansion or a finite-cylinder warm start:

```sh
julia --project=. test/infinite_block_initializer.jl
julia --project=. test/infinite_number_background.jl
julia --project=. scripts/infinite_random_cell.jl armchair 2 16 64
```

This example fixes mean physical density to 16/36 = 4/9. A normalized library
`random_mps` supplies an entangled finite cell; deliberate unit-dimensional
boundary bonds repeat independent cells before official VUMPS expansion.
Every physical spin remains an individual MPS site. Shifted/scaled site QNs
encode the filling and leave physical spin matrices unchanged. Zero virtual
QN flux therefore does not imply half filling. Output records the physical
number background explicitly and uses `rho4of9` in checkpoint names. These
are alternative variational candidates, not a certificate of the globally
preferred filling or a two-dimensional phase.

For an original QN infinite-state payload, transfer eigenvectors can also be
classified by their virtual physical-number difference:

```sh
julia --project=. test/transfer_charge_labels.jl
julia --project=. scripts/infinite_transfer_charge_audit.jl RESULT.jls
```

The dense eigenvector basis is matched to the original site and bond QN labels.
Degenerate vectors may mix charge blocks, so the output saves all weights.
This identifies U1 transfer channels, not microscopic CGS charges or a gap.


Infinite unit cells now support an explicit `cell_slices` parameter; the
historical default remains two. The checked three- and six-slice definitions
are saved in `geometry/infinite_*_cell{3,6}_*.json`, regenerated with
`julia --project=. scripts/export_infinite_periods.jl`, followed by
`python scripts/expand_infinite_site_records.py` and
`python scripts/strategy_geometry_audit.py`. They use individual
physical spin sites and the same exact endpoint exchange coefficients.
For example, a separate six-slice density-4/9 armchair width-two candidate is
`julia --project=. scripts/infinite_random_cell.jl armchair 2 48 64 6 parallel`.
This command starts a new branch; existing two-slice checkpoints are preserved.
Local solver residuals are audited, and larger cells have unique tags and
checkpoint names. Mean slice entropy is descriptive: inspect every cut and
modulation before fitting. Transfer lengths in slices equal `cell_slices`
times the cell transfer length; energy per vertex uses `2*w*cell_slices`.
These alternatives have not yet supplied a converged 2D phase diagnosis.

Infinite allocated bond dimensions include the unit-cell wrap bond in both AL
and AR. The pinned finite-library `maxlinkdim` omits that bond; historical raw
counts are retained alongside `results/infinite_bond_dimension_audit.json`.
Reaudit saved tensors with `julia --project=. scripts/audit_infinite_bond_dimensions.jl`.
Requested cap, allocated dimension and effective Schmidt rank are different.
Full slice-cut probabilities are saved by new measurements. A raw energy above
a trial state supports exclusion only after canonical/normalization checks;
fixed-density comparisons also require a trial with compatible filling.

The validated cell initializer accepts an optional site ordering after the
update algorithm. For example,
`julia --project=. scripts/infinite_random_cell.jl armchair 2 26 32 3 parallel matter_first`
tests physical density 13/27 in 54 individual spins. This is a separate candidate,
not a preferred-filling assumption. At armchair width two, matter-first ordering
reduces the maximum operator span from 31 to 17 sites; its internal entanglement
requirements still need comparison. Both orderings' nonlocal center operators
pass independent microscopic complex-field checks. Original star-order runs
and their checkpoints are preserved. The conditional ordinary D(Z3) filling
constraint and initializer accessibility are in `results/primitive_filling_audit.json`.

The scientific re-audit also checks convergence and admissible fillings without
restarting a production optimization:

```sh
OPENBLAS_NUM_THREADS=1 python scripts/primitive_filling_audit.py
OPENBLAS_NUM_THREADS=1 python scripts/number_bound_audit.py
OPENBLAS_NUM_THREADS=1 python scripts/analyze_infinite.py
OPENBLAS_NUM_THREADS=1 python scripts/convergence_gate_audit.py
```

The number bound is an exact hard-core-spin operator inequality, not a free-boson
approximation. Combined with a valid 2D star-product trial, it rules out extreme
bulk densities; it does not select a density or establish topological order.
The convergence report treats missing evidence as unknown and keeps recorded
optimization criteria separate from global-sector and thermodynamic claims.
Odd-width half-filled rows have fractional number per primitive helical
translation, so parity and symmetry branches must be checked before pooling
circumferences. See `results/scientific_strategy_reaudit.md` for the current
interpretation of all preserved results. Future large jobs should run serially
while existing workers cause substantial memory compression.

Matter-triplet measurements and their historical interpretations are retained
in `results/archive/triplet_provenance/` and the original raw result files.
They are excluded from core diagnostics, phase interpretation, convergence
criteria, and entropy fits. No distinguished microscopic role has been derived.
