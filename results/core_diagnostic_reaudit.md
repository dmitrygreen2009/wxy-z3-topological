# Core diagnostic policy re-audit

The matter S+S+S+ one-point function and paired matter-triplet correlations are
removed from the core scientific pipeline. CGS invariance is an algebraic
property; it does not by itself establish a distinguished physical diagnostic
for this model. No such derivation has been supplied here.

Superseded: the interpretation that a measured triplet one-point function or
disconnected plateau provides relevant U(1)-breaking evidence for the phase
study; comparisons using triplet correlation amplitudes to judge transfer-mode
relevance; recommendations to extend triplet measurements as a primary test.
No convergence or entropy-fit acceptance criterion may depend on those data.

Preserved: original raw triplet JSON measurements and tests remain available
for provenance. Pre-removal versions of the affected audits and summary are
in `results/archive/triplet_provenance/`. Matter blocks in exact CGS operators
remain mandatory and are unrelated to the discarded observable. Benchmark
spectra, Hamiltonian, geometries, energy optimizations, full-state Schmidt
spectra, conserved-number correlations, and transfer eigenvalues are unchanged.
`results/torus_hard_gate_audit.md` remains untouched.

The neutral number-selection rule for winding operators and the measured
transfer eigenvector number labels remain independently supported. Their use
as geometry-specific flux evidence still requires the translated-loop audit.
No accepted topological entropy fit exists. No 2D phase is established or
excluded by the current data.

Future finite/infinite core observable drivers and the microscopic-observable
summary omit triplet measurements. Permanent `test/core_diagnostic_policy.py`
checks those exclusions and archival preservation. Existing running drivers
may have loaded older code; any subsequent legacy triplet output is archival
only and cannot enter interpretation or convergence decisions.
