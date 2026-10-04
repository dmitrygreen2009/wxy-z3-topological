# Scientific strategy re-audit

This audit uses the existing exact-Hamiltonian source, saved spectra, MPS results,
geometry records, and convergence logs. Validated ED/DMRG/VUMPS calculations were
not restarted. New static embedding/homology checks and local symmetry identities
were added; ongoing runs and all valid physical checkpoints were preserved.

**Scientific status: the current numerics neither establish nor exclude a
2D Z3 topologically ordered phase. They also do not establish 2D gaplessness
or spontaneous physical U(1) order. No current entropy intercept is an accepted
thermodynamic topological gamma.** This audit is the current interpretation of
the historical results and takes precedence over shorthand in older progress notes.

## Hamiltonian and geometry

The microscopic implementation has the specified complex DFT3 W/sqrt(3), J=1,
local three-spin matter triplets, and one spin per physical gauge edge. Raising
and lowering amplitudes are conjugate partners. The corrected torus multiplicities
4,1,1,4 remain validated; the historical hard-gate audit is unchanged. Neither
fourfold degeneracy on that tiny quotient nor its finite-sector excitation gap
identifies the thermodynamic topological order.

The physical embedding is a1=(sqrt(3)/2,3/2), a2=(-sqrt(3)/2,3/2), B-A=(0,1).
Zigzag identifies (x,y) with (x,y+w), has circumference sqrt(3)*w and axial
coordinate t=x. Armchair identifies (x,y) with (x+w,y+w), has circumference
3*w and axial coordinate t=x-y. Axial distance per slice is 3/2 and sqrt(3)/2,
respectively. Both endpoints use the explicitly recorded local leg. Parallel
edges remain distinct. Open ends retain dangling physical gauges: Nspin is
9*L*w+w (zigzag) or 9*L*w+2*w (armchair), rather than the periodic bulk count.

The independent static audit passed **41 finite embedded cuts and 16 infinite
incidence/cut tables**. It checks physical positions rather than merely agreement
of two code orderings. In units of the axial step, zigzag matter A/B positions
are t and t+1/3; gauge legs 1,2,3 are t+1/6,t-1/3,t+1/6. Armchair matter
positions are t; gauge positions are t,t-1/2,t+1/2. The half-open cut is at
L/2. Integer-scaled coordinates confirm both even and odd recorded cuts and
both infinite slice cuts. Matter triplets remain wholly on one side and every
physical spin belongs to one side. The infinite cell contains two slices,
4*w vertices and 18*w physical spins. Saved hopping tables independently match
the unfolded honeycomb incidence, including off-cell exchanges.

The new winding audit exports 221 explicit circumference cycles. Importantly,
the zigzag width-one parallel two-edge cycle is **noncontractible in the original
honeycomb embedding**: its lift ends at displacement (0,1). Calling it a local
two-cycle describes its support on the narrow quotient, not contractibility.
The extra local S3 reflection is real and independently verified; it can produce
sector degeneracies peculiar to this collapsed geometry. Width one cannot supply
reliable 2D intercept or gap evidence. Wider hexagon charges and circumference
cycles must be distinguished by their unwrapped winding.

## Filling and actual symmetry sectors

The benchmark number blocks are correct. They do not fix the ground-state
number block on open cylinders. Complete small-cylinder ED establishes zigzag
L2,w1 minima at N_up=9,10 and armchair minima at 8,12, rather than half-filled 10.
The independent block spectra and off-half ITensor comparison remain valid.
Their gaps are small-system, fixed-sector gaps.

Larger number-sector searches are variational, not complete ground certificates.
For example zigzag L4,w1 has a saved N_up=17 candidate at -15.0337476103,
below the particular half-filled/trivial-cycle candidate near -15.0328573215.
Armchair L4,w1 likewise has off-half candidates below the current neutral
variational state. These comparisons do not lower-bound the unexplored neutral
sector or prove a bulk density shift. A boundary shift of O(w) and an extensive
shift of O(L*w) must be distinguished with length and density profiles. Exact
particle-hole conjugation guarantees equal mirror-sector energies; discrepancies
between independently optimized partners are solver effects, not physical
particle-hole splitting. The half-filled infinite QN initialization is an
assumption about density, not a demonstration of optimal bulk filling.
For even finite L, zigzag has Nspin=(9L+1)*w: odd widths have an odd
physical spin count and a nearest-half sector with total Sz=+/-1/2, whereas
even widths can have Sz=0. Both must be tracked; parity-dependent boundary
or sector corrections cannot silently be absorbed into a common intercept.
The two-slice infinite cells have even Nspin=18*w and imposed Sz=0.

A concrete sector-control failure was found without invalidating its numerical
energy. The L4,w1 seed7214 run passed its initial q=1 projection checks but
finished with <U_t> nearly 1 on every cycle. Its final energy -15.0328573116575
therefore does **not** measure a nontrivial-sector ground energy. Two-site
optimization/truncation does not enforce the multi-site CGS charge merely because
it was initially projected, even with zero noise. The original result and loop
measurements are retained and annotated. All future projected-initializer results
must distinguish requested initial charge from independently measured final charge.
There is no reason to rerun its validated energy as an unconstrained variational
result. A conditional-sector calculation would require additional enforcement
and repeated purity checks, not reusing its initial label as a final constraint.

Exact microscopic CGS charges are not automatically identified topological
anyon/MES labels. All nine local affine exponent transformations have now passed
independent 64-ket commutator, order-three, and unitarity checks. A winding
operator can be a useful candidate flux diagnostic, but its logical meaning and
which basis minimizes the entropy still need identification. The narrow S3
charges particularly cannot be treated as thermodynamic anyon counting.

## Convergence and correlation lengths

Stationarity, state accuracy, sector-ground optimization, and thermodynamic
convergence are separate questions. Small energy variance and small sweep drift
can describe an excited sector, a trapped state, or an arbitrary degenerate
superposition. Actual discarded weight is controlled by the achieved bond dimension;
a requested cutoff of 1e-11 does not imply discarded weight is below that cutoff
when maxdim binds. Width-three chi512 logs still have truncation errors around
1e-5 to 1e-4, above the recorded 1e-8 target. Several narrow conditional states
are stable in energy/entropy, but this does not certify a global cylinder minimum.

The stored criteria remain: sweep energy drift 1e-8, energy-per-vertex chi drift
1e-6, chi entropy drift 1e-3, refinement entropy drift 1e-4, discarded weight
1e-8, and length entropy drift 1e-3. A fixed sweep count is a work stage, not a
convergence certificate. Criteria not checked or lacking a compatible comparison
must be marked unknown, not passed. No tolerance was relaxed in this audit.

The analytic entropy and transfer tests validate the numerical pipeline. Transfer
lengths use the library transfer matrix for a real infinite state, include all
charge sectors through measurement-only densification, and retain complex leading
eigenvalues and residuals. One transfer cell has two axial slices. Thus xi_cells
and xi_slices differ by two; physical axial xi is 3*xi_cells (zigzag) or
sqrt(3)*xi_cells (armchair). This is an axial correlation length, not an independent
measurement of the transverse 2D correlation length or a spectral gap.

Independent canonicalization defines a consistent state from the stored AL
representation. It does not optimize the state or certify that the original,
slightly inconsistent AL/AR/C triples were identical. Energy/entropy changes and
the inherited VUMPS residual must remain visible. The width-one QN zigzag chi128
fixed-point measurement gives mean entropy 1.27225627787 and xi_slices 17.0616437381,
with variational residual 3.2566e-5. Its increased xi from chi32 and chi64 is not
proof of a 2D critical mode: its cycle charges are mixed and the geometry has the
extra local S3 degeneracy. Conversely, finite xi at finite chi in a product-like
or unconverged state is not evidence of a gapped physical ground state.

Two random dominant-transfer starts detect some repeated fixed points; the
reported detected rank is not a certified full fixed-point multiplicity. A unit
modulus peripheral eigenvalue can reflect a cat, multiple sectors, or a larger
translation period. Physical connected correlations within an identified state
branch, gauge-invariant operator overlaps, and chi convergence are needed before
turning transfer eigenvalues into a bulk conclusion. Different initializer tags
are now kept separate in chi convergence comparisons.

An exact trial-state bound exposes particularly misleading small residuals.
The disjoint-A-star product gives E0_cell <= -4.8061842*w, conservatively rounded;
independent 18-spin expectation error is below 6e-15. Armchair width-two cap2 and
cap4 points violate this bound despite tiny reported residuals and are excluded
as ground-state candidates. Their raw records are retained.

## What entropy fits can and cannot establish

All primary entropies include matter and gauges in the same pure state. Gauge-only
mixed-state constants are not used. Nevertheless, full-state entropy can still
carry state-selection, global-number, local-degeneracy, symmetry-restoration,
and end-to-end contributions.

For a sufficiently long gapped cylinder in an identified minimally entangled
branch, the asymptotic single-cut constant is ln(D/d_a). For Abelian D(Z3),
D=3 and d_a=1, hence ln(3). Coherent sector superpositions can add entropy, so
a finite-system exact energy eigenstate or an eigenstate of a particular loop
need not be the MES appropriate for this test. Lowest finite energy alone is
therefore not a sufficient entropy-input selection rule. The MES argument and
aspect-ratio issues are discussed in [Jiang, Wang and Balents](https://arxiv.org/abs/1205.4289),
while [Zhang et al.](https://arxiv.org/abs/1111.2342) demonstrate state dependence
and the need for additional data to identify anyon statistics.

The same nominal L across increasing w does not hold the physical aspect ratio
fixed: nominal Lx/Ly is sqrt(3)*L/(2*w) for zigzag and L/(2*sqrt(3)*w) for armchair.
At L4,w4 these are about 0.866 and 0.289. Neither family has established
end-insensitive bulk entropy at those widths. Finite and infinite points with
different filling, actual CGS content, or initializer branches must not be
combined as if they were the same thermodynamic state. Nor should axial xi from
one mixed width-one branch certify length convergence of a different finite q=0
branch or transverse convergence at another width.

U(1) number conservation is exact. A zero one-point charge-three matter expectation
in a QN ansatz is enforced by that symmetry and cannot exclude spontaneous order.
Raw <S+> can also be gauge variant; it is a variational/coherence diagnostic,
not a physical order certificate. Gauge-invariant triplet correlations and their
length/width scaling are needed. Continuous symmetry breaking can produce
logarithmic entropy terms rather than the assumed constant-only correction;
see [Metlitski and Grover](https://arxiv.org/abs/1112.5166). This is a competing
possibility to test, not a claim that this model breaks U(1).

Existing OLS errors measure scatter about a chosen finite-data line. They do not
include uncontrolled chi, length, sector, and circumference systematics. Negative
or unstable apparent gamma does not exclude Z3 order, and a coincidental intercept
near ln(3) would not establish it. Descriptive fits are retained separately;
accepted topological fits are currently empty. Even a reliable ln(3) constant
would strongly constrain total quantum dimension but would not alone distinguish
untwisted D(Z3) from every other order with the same D.

## Revised continuation priorities

1. Preserve ongoing exact-Hamiltonian runs and resume their newest valid states.
   Increase dimensions where actual truncation/entropy drift requires it; do not
   restart a converged benchmark or treat a completed stage as convergence.
2. Finish neutral/off-neutral candidate refinements and wider-cylinder number
   searches. Add bulk/end density profiles to distinguish boundary charge from
   a competing bulk density. Preserve both sides of particle-hole comparisons.
3. Measure actual contractible and noncontractible CGS content of existing
   wider-cylinder states before pooling entropy data. The failed pure-q1
   initialization is an algorithm audit, not a sector-energy result. Extra effort
   on the width-one quotient must serve debugging, not dominate the 2D inference.
4. Compare independent infinite initializations at width >=2, retain their branch
   identity, and verify physical connected correlations. The fixed two-slice
   infinite unit cell can miss larger-period phases; finite bulk profiles and,
   if indicated, a validated larger cell must test that restriction.
5. Establish compatible finite-length and chi plateaus for each usable width.
   Then enlarge circumference, exclude the special narrow quotient, and compare
   separate wrapping fits with smallest-width removal and systematic stability.
   Examine MES/flux basis selection and any quasi-degenerate branches explicitly.
6. Only after these conditions hold may a stable gamma and finite physical xi be
   used as positive evidence. If available resources cannot reach those conditions,
   the scientifically justified endpoint is an audited inconclusive result with
   explicit achieved limits, not an unstable-intercept verdict.

Concrete fixes from this audit are branch-separated analysis, explicit fit
ineligibility, source-checksum correction for a charged input alias, actual-charge
annotation of the migrated initializer, and machine-readable winding definitions.
Legacy warm-fit caches lack complete input fingerprints; their optimized states
remain valid variational states, but exact initializer provenance is incomplete.
This is a reproducibility gap to correct for future caches, not evidence that the
validated Hamiltonian energies should be discarded or recomputed.

The new width-two circumference measurements show different nearly pure charges at different axial rings: zigzag probabilities are predominantly q=2 at t0 and q=0 at t1–t3; armchair is predominantly q=2 at t0,t1 and q=0 at t2. Operator-application norm errors are below 1e-12. This is useful evidence about saved-state and boundary/flux structure, not proof of a homogeneous topological MES. The small residual charge admixtures and bond/length dependence still require convergence. Density profiles are saved alongside these new observables.

## Follow-up audit of analysis coverage and physical correlations

The finite analysis readers had omitted completed number-scan candidates whose
metadata were stored inside each record, rather than at the outer file level.
This biased the displayed lowest-energy envelope toward the reference filling.
A shared reader now includes those candidates, PH-transformed states, and exact
preserved result snapshots; a permanent fixture regression checks the omission.
No production state or validated benchmark was recomputed. Entropy fits remain
descriptive and ineligible for topological inference. The plotting code also
separates initializer and number-sector branches instead of joining unrelated
points into apparent convergence curves.

For the consistently canonical width-one QN cap128 state, the leading full
transfer magnitudes are 1, 0.8893877 (twice), and 0.7396177. New winding-loop
connected contractions at 8 and 16 cells give a descriptive per-cell decay
ratio about 0.7384, roughly following the fourth eigenvalue rather than the
leading subdominant pair. Charge-three matter contractions are much smaller.
Neither operator measurement identifies the leading eigenvectors' quantum
numbers, and neither certifies a physical gap or ground-state order. Raw
contractions and a separate comparison are saved in
`results/transfer_operator_coupling_audit.json`. Existing validated one-point
values were reused; legacy records lack an original payload fingerprint, which
is distinguished from the checksummed payload and one-point JSON used by these
new contractions. Future one-point measurements fingerprint their payloads.

The filling audit also tests an independent fixed-density infinite initializer.
Ten exact-state checks validate repeating an entangled finite cell with an
explicit one-dimensional boundary bond and official library canonicalization.
Twenty-one checks validate the physical density encoded by shifted/scaled site
QNs, including importing both KrylovKit and ITensors and verifying that `Sz`,
`S+`, and `S-` retain their physical spin matrices. A zero shifted virtual flux
is not a statement of physical half filling. New output records the background
explicitly; legacy dense measurement tensors without that provenance remain
unclassified rather than being assigned a density retrospectively. The new
armchair width-two density-4/9 run is an independent candidate, motivated by
the finite number scans, and provides no preferred-filling certificate yet.

The first production expansion exposed a boundary-arrow incompatibility not
covered by the ten state-only initializer checks. Both initial failed logs are
preserved. The repair changes only the direction of a one-dimensional zero-QN
right-canonical boundary index, consistently in AR and C, using the library's
index operations; it changes no amplitudes or physical spin operators. Thirteen
checks now include successful official subspace expansion across that boundary.
The failed run produced no optimized scientific state and was retried from its
same seeded initializer after this concrete interface correction.

New consistent-state measurements of the completed armchair width-two QN
cap16 stage give energy -9.525782017201935, canonical error 2.04e-15 and
spatial entropy below 1e-30. Recanonicalization changes its energy by only
3.2e-14. The compatible isolated-A-star trial energy is at most -9.6123684:
this particular state is excluded as a sector ground state, despite its
inherited 3.15e-6 solver residual. Its zero transfer length is a slice-product
artifact, not a phase conclusion. No optimization was repeated for this audit;
the newer QN stage continues. See `results/armchair_qn_initializer_trapping_audit.json`.

A new eigenvector measurement using the original width-one QN payload resolves
part of the transfer interpretation. The two leading subdominant eigenvectors
at magnitude 0.8893877 have virtual physical-number differences -1 and +1,
respectively, with flux-weight purity within rounding of one and eigenvector
residuals below 9e-14. The next mode at 0.7396177 is neutral and has transfer
length about 6.63 slices, compatible with the descriptive winding-correlation
decay about 6.60 slices. Thus the full length 17.06 does not control the neutral
winding or charge-three matter operators in this U1 ansatz. The classification
does not identify microscopic CGS charges or turn this unconverged width-one
state into evidence for a two-dimensional gap or phase. Ten analytic matrix
checks validate the QN-basis flux-weight bookkeeping; all transfer contractions
and eigensolves remain library operations. Original optimizations were reused.

New saved-state armchair L4,width2 contractions show that the Nup=31,33 and
38 candidates also differ in approximate circumference-charge patterns: their
most-weighted charges are (0,1,1), (1,2,0) and (2,2,0), respectively. Purity
errors remain about 1e-4--1e-3, so these are not exact sector certificates.
Relative to Nup=38, the central half-open interval [1,3) contains about -5.034
of the total -7 number change at Nup=31, and -3.515 of the total -5 at Nup=33.
These short, unconverged states do not distinguish bulk filling from end and
sector effects. See `results/armchair_L4_w2_number_and_cycle_audit.json`; the
raw circumference probabilities and density profiles are preserved separately.


A primitive translational cell contains two stars, six matter spins and three
owned shared gauge spins: nine physical spins. The explicit record is
`geometry/primitive_honeycomb_cell.json`. Thus physical half filling has
U1 number 9/2 per primitive cell, whereas densities 4/9 and 5/9 have integer
numbers 4 and 5. This is a substantive constraint on interpreting a fixed
half-filled calculation, not a reason to change its measured results.

For a gapped phase preserving U1 and primitive translations, the usual
background-anyon formulation requires fractional anyon charge equal to the
filling modulo one. See Cheng et al.,
[Translational symmetry and microscopic constraints](https://arxiv.org/html/1511.02263v3),
section VI.2, and Zaletel and Vishwanath,
[Constraints on topological order](https://arxiv.org/abs/1410.2894).
Applying this formulation to ordinary D(Z3) is our inference: all its anyons
have fusion order dividing three, so charge additivity permits only 0, 1/3
or 2/3 modulo one, never 1/2. Under the formulation's symmetry assumptions,
a symmetry-preserving gapped half-filled D(Z3) candidate is incompatible with
this constraint. This conditional argument does not exclude the model's
realizing Z3 order at another filling, with broken translations/U1, or with
additional order. Exotic translation actions permuting anyon types have not
been independently audited here; no blanket exclusion is asserted.

The two-slice infinite ansatz is also a restriction. Translation can permute
minimally entangled sectors without breaking local bulk translation symmetry;
a three-sector cycle may require a three-slice representation. A six-slice
cell can accommodate both periods two and three. Current two-slice fixed points
therefore do not exhaust possible minimally entangled states, and variation
among slice observables must be tested before inferring physical symmetry
breaking. Larger-cell implementation and comparisons remain required.

Local eigensolver tolerances were previously requested without saving all
achieved local eigenpair residuals. New VUMPS logging saves ConvergenceInfo
and independently evaluates norm(M*v-lambda*v)/norm(v). The stopping condition
now also requires every local solve to meet its unchanged requested tolerance.
Six analytic solver checks and five comparisons against the official parallel
one-iteration driver pass; the latter reproduces the same tensors and records
36 local solves. Existing independent benchmark validations remain valid.
The official simultaneous-update VUMPS branch is additionally checked on an
exact two-site ferromagnetic Ising state (six checks). Its environments are
reused within an iteration; it is not a new tensor-network implementation.
An initial onsite-only probe exposed an unsupported range-one MPO interface;
its failed log is preserved and is not scientific evidence. Existing WXY
calculations contain the supported two-site exchange terms. Active original
jobs continue unchanged; new continuations use the stronger achieved-residual
audit. A converged local solve still does not certify the global ground state.


Periods three and six now have explicit microscopic pair tables in the
repository, with the historical two-slice default preserved. Independent
embedded-coordinate checks pass for all 40 infinite tables. The generic
factory passes 216 spin-incidence checks, 12 coherent-spin exchange-energy
checks and nine analytic transfer-period tests. Seven end-to-end six-slice
product-state measurement checks verify all six spatial cuts, energy per
vertex, trial-bound normalization, transfer-length units and checkpoint names.
The known product fixture is excluded from scientific fit input. Resume and
recanonicalization use the saved state cell size, and number-background and
initializer branches remain separate in analysis. Production calculations
with larger periods remain necessary; implementation validation is not
variational convergence or evidence of topological order.


New central retained-basis residuals for the unchanged armchair L4,width2,
Nup31,chi256 state are 0.0040861, 0.0031193 and 0.0030396 on bonds 39--41.
Each projected Rayleigh energy matches its saved -30.6119657548134 energy
within 6e-14. Five exact two-spin tests validate the library ProjMPO diagnostic.
These are saved-state stationarity residuals, including truncation and basis
limitations; they are not the unrecorded pre-truncation Krylov residuals.
They provide concrete additional grounds for continuing this underconverged
candidate, without changing the validity of its energy as a variational bound.
The continuation interface now records configurable Krylov iteration controls;
the historical defaults remain unchanged and requested tolerances are not
loosened. Ground-state and entropy convergence still require separate tests.


The infinite cell translator advances primitive a1 by the recorded number of
slices. It generally also shifts around the circumference; it is not always
a pure perpendicular axial displacement. New geometry metadata records both
components explicitly. The saved physical axial transfer length correctly
uses the perpendicular projection (1.5 per zigzag slice, sqrt(3)/2 per armchair
slice), whereas a complex transfer phase can also contain transverse momentum.
Fixed-site correlations across cells follow this helical displacement. This
clarifies the convention without changing existing numerical transfer lengths
or treating them as a two-dimensional spectral-gap certificate.


An additional interpretation correction separates a raw infinite canonical
expectation from a normalized consistent-state variational bound. At sizeable
AL*C versus C*AR mismatch, a center-based expectation cannot by itself prove
that a physical ground candidate lies above a trial state. Future measurement
and saved-data analysis explicitly check canonicality, both isometries, center
normalization and transfer normalization at 1e-10 before supporting that
exclusion. A fixed-density comparison additionally requires the isolated-A-star
trial to have a compatible filling (the construction covers densities 1/3
through 2/3). Raw energy comparisons are preserved. The previously audited
armchair QN cap16 exclusion remains supported: the consistently recanonicalized
state has error 2e-15, changes energy by 3e-14, has compatible half filling,
and lies about 0.0866 above the explicit trial. None of these checks certifies
the global ground state or a thermodynamic entropy intercept.


The same conditional filling argument allows density k/27, not only integer
number per primitive cell. The closest permitted densities to half filling
are 13/27 and 14/27, with fractional primitive-cell charges 1/3 and 2/3.
Thus the density-4/9 candidate does not exhaust the symmetry-compatible
alternatives. The current integral-number cell-product initializer cannot
initialize 13/27 in a width-two two-slice cell; three or six slices permit
exact integral cell numbers 26/54 or 52/108. This is an initializer limitation,
not a claim that every possible noninjective two-slice MPS excludes that mean
filling. The arithmetic and accessibility table are stored in
`results/primitive_filling_audit.json`. Theory-compatible densities are test
candidates, not presumed global ground fillings.


The dense transfer audit exposed an allocated-bond reporting issue. The pinned
finite-library `maxlinkdim` scans n-1 bonds and omits the infinite unit-cell
wrap bond. The completed zigzag width-two cap32 parallel state was reported
as chi16 but its wrap bond has dimension29, giving a full virtual transfer
space of dimension841. Every eigenvalue of that space was independently
obtained through official TransferMatrix actions and dense diagonalization:
there is one peripheral eigenvalue at the stated 1e-9 tolerance, the maximum
residual is 1.74e-13, and xi=0.117676310549364 cells agrees with scalar Krylov.
This checks the spectrum of the same AL state, not a physical bulk gap.

Future logging, filenames and memory estimation inspect all periodic AL and AR
links. Eight index-level checks validate the correction, including a larger
wrap bond and a finite-chain control. A saved-state audit corrects allocated
bond counts without repeating optimizations or replacing historical raw data.
Allocated dimension, requested cap, and effective Schmidt rank are distinct.
Future measurements also retain full probabilities at every spatial slice
cut. Near-zero Schmidt support can leave unused auxiliary transfer modes;
operator coupling and boundary weights matter when interpreting those modes.
The already documented two-start peripheral detection remains a lower-bound
method except where a complete small-space dense audit is saved.

Independent complex-product checks validate the infinite Hamiltonian against
explicit vertex/endpoint records (12 checks), rather than only its coherent
sum of W entries. The official nonlocal VUMPS center operator also matches
independent microscopic local 2x2 fields for both endpoint families and both
site orderings (48 checks). This validates complex coefficients, Hermiticity,
periodic terms and local effective operators without rerunning benchmarks.
The library pads all term MPOs to their common maximum range, satisfying its
constant-range assumption. No Hamiltonian approximation was introduced.

The completed zigzag width-two parallel state was remeasured from the same AL
tensors with consistent centers, without optimization. Canonical error fell
from 4.50e-5 to 6.83e-12; energy changed by 9.29e-7 per two-slice cell and
entropy by 3.05e-6. Its full dense 841-dimensional transfer spectrum agrees
with the Krylov correlation length (0.11767631055 cells), with eigenpair
residuals below 1.74e-13 and one peripheral eigenvalue at tolerance 1e-9.
This confirms this measurement pipeline, not convergence of the ground state
or a physical 2D gap. The complete finite-dimensional transfer builder also
passed known primitive-channel and GHZ tests, including exact multiplicity.

Additional density branches remain planned rather than launched during heavy
compressed-memory pressure. Existing expensive calculations continue, and
smaller completed results and their local checkpoints remain preserved.

Checkpoint-resume source inspection found a future reproducibility hazard:
finite resumes unconditionally used Krylov dimension 12 and default one
local iteration, even for a saved continuation specifying dimension 32 and
ten iterations. The resume driver now retains recorded controls, including
the refinement-specific dimension, for both the resumed sweep and subsequent
bond stages. No active run was interrupted or restarted. Infinite resumes
also explicitly retain seeds and record their source checkpoint provenance.

The density strategy must not prioritize only values closest to one half.
The short armchair charged-state profiles differ substantially between
central rows and ends, while their energy splittings are comparable to
remaining optimization errors. Neither their total filling nor the densities
13/27 and 14/27 certify a bulk optimum. The complete conditional k/27 grid
is recorded; broader density searches, length dependence, and particle-hole
partners remain required wherever competitive energies emerge. A gapped
plateau requires stability against number changes, not merely one fixed-N
optimization. No density or symmetry branch currently certifies the 2D phase.
