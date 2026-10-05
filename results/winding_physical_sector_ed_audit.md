# Sector-resolved ED at the audited small-cylinder fillings

This is a new check needed for direct winding-sector optimization, not a
rerun of the validated torus benchmark or existing unrestricted spectra.
The independently derived full eight-state real exchange form is assembled
in each exact number-and-winding block using the recorded shared-edge and
endpoint leg definitions. No matter state or shared gauge spin is discarded.
First, all six N2 sector ground energies agree with the earlier original-basis
projector ED within 1e-9. The new higher-filling results are:

| Family | Spins | N_up | S_z | Filling | q | Block dimension | Ground energy |
|---|---:|---:|---:|---:|---:|---:|---:|
| zigzag | 19 | 9 | -0.5 | 9/19 | 0 | 31460 | -7.696107978815055 |
| zigzag | 19 | 9 | -0.5 | 9/19 | 1 | 30459 | -7.651758391631642 |
| zigzag | 19 | 9 | -0.5 | 9/19 | 2 | 30459 | -7.651758391631634 |
| armchair | 20 | 8 | -2 | 2/5 | 0 | 42000 | -8.037741645057547 |
| armchair | 20 | 8 | -2 | 2/5 | 1 | 41985 | -7.495281096632541 |
| armchair | 20 | 8 | -2 | 2/5 | 2 | 41985 | -7.495281096632553 |

These numbers are fixed-sector minima at previously independently audited
finite ground fillings. They do not select a thermodynamic filling. The three
sector dimensions sum to the complete fixed-number dimensions 92378 and
125970. The minimum across winding sectors agrees with the existing,
independent original-basis block ED to tolerance 1e-9; that comparison uses
preserved results and does not rediagonalize the original Hamiltonian.

The solver is real sparse SciPy ARPACK, seed 7258, k=8, ncv=80, tolerance
1e-12 and maxiter=10000. Every returned vector is independently orthogonalized
within its degenerate cluster and checked against the complete block H.
The maximum eigenvector residual is 5.3e-12. Hermiticity is exact in the
assembled real sparse representation. Energies, per-vector residuals, Gram
errors, runtimes, explicit winding weights and source provenance are saved
in `winding_physical_sector_ed.json`; fixture cross-checks are saved separately
in `winding_real_sparse_fixture_ed.json`.

Scalar ARPACK does not certify complete degeneracy multiplicities. No claim
of complete degeneracy is made from its returned eight eigenpairs. The
earlier complete block-ED audit establishes that the q0 global ground state
at each of these two fixed fillings is unique. In the armchair nontrivial
blocks, four orthogonal ground vectors are returned; a chosen optimized MPS
must still establish entropy stability and any claimed MES separately.

The finite winding-sector gaps are approximately 0.0443495872 (zigzag) and
0.5424605484 (armchair). These are exact small-cylinder symmetry-block energy
differences, not deconfinement energies, asymptotic topological splittings,
or a determination of the 2D phase. No entropy from this audit is fitted.
