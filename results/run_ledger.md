# Run ledger

All energies and entropies in JSON files are measured results, not estimates.
The final achieved caps are in the records; a filename's chi is the planned
cap, which an interrupted run may not have reached.

* `ed.log`, `ed.json`: independent SciPy fixed-number ED, passed all ground
  energy tests and resolved the tiny torus's additional degeneracies.
* `validation.log`: first ITensor benchmarks (ITensorMPS 0.4.2), passed.
* `validation_pinned.log`: rerun after adding the compatible infinite package
  selected ITensorMPS 0.3.45, passed. This is the committed environment.
* `tests_final.log`: lattice, spatial-cut ordering invariance, infinite
  physical-incidence, and MPO Hermiticity regression checks, passed.
* `batch.log`: original axial ordering, L=4, w=1–3 for both families,
  planned and achieved chi=128. These are exploratory, unconverged wider data.
* `refined_batch.log`: original axial ordering, L=6. Zigzag w=1 reached 256;
  zigzag w=2 was interrupted during the 256 stage after star ordering reached
  a lower energy at 128 much faster. Its JSON/checkpoint contains completed
  lower-cap stages. The interrupted stage remains in the log.
* `extended_batch.log`: original axial zigzag L=10, w=1 reached a completed
  256 stage and two logged 512 sweeps. Interrupted to replace the expensive
  axial ordering with star ordering. No 512-stage convergence is claimed.
* `infinite_zigzag_w1.log`: initial product-state VUMPS attempt, interrupted
  after low-dimension stagnation. Not used for correlation lengths.
* `infinite_zigzag_w1_v2.log`, `infinite_armchair_w1.log`: rotated product-state
  VUMPS, spatial cut matched to the finite convention, planned chi=32.
* `ordered_batch.log`: improved star ordering, L=6, w=2 for both families,
  planned chi=128, with noiseless refinement.
* `zigzag_L6_w2_chi256_star.log`: checkpoint continuation with four Julia
  threads for ITensor block-sparse operations; BLAS remains single threaded.
* `scaling_batch.log`: improved star ordering, L=4, w=1–3, both families,
  planned chi=256, four Julia threads.

The two orderings differ only inside spatial slices and have identical sets
of physical sites to the left of the measured spatial cut. No circumference
supersite, clock approximation, or matter trace was used.
