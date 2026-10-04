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
