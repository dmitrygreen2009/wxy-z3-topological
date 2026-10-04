# Scientific and repository workflow

Use the exact microscopic WXY spin Hamiltonian. Each shared gauge edge is one
physical spin; keep the local matter spins in full-state entanglement measurements.
Use ITensors/ITensorMPS for production tensor networks and independent ED where
practical. Preserve results/torus_hard_gate_audit.md as the historical audit.
The corrected full N_up=9 torus spectrum has multiplicities 4, 1, 1, 4 for
the first ten levels. Validate multiplicities and eigenvector residuals, not
just the ground energy.

After each major validated milestone, commit source code, regression tests,
lattice definitions, environment requirements, analysis scripts, documentation,
audit notes, and small reproducibility results. Push to the current GitHub
branch. Do this at a natural checkpoint without interrupting a running
calculation. Keep large MPS/checkpoint files, large generated binaries,
package caches, and temporary scratch files out of Git; update .gitignore.

If pushing fails, preserve every local commit and working file, record the
exact Git error, and continue scientific calculations where possible. Before
the overall project finishes, verify that GitHub contains the complete
reproducibility workflow and final numerical summaries.
