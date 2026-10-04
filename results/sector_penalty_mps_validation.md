# Small-cylinder exact sector-optimizer validation

The exact commuting projector penalty produced converged q=0 and q=1 states for both cylinder families in independently ED-checkable N_up=2 blocks. This validates sector isolation and the optimizer. These low-filling fixtures are excluded from ground-phase summaries and entropy fits.

| Family | Spins | N_up | S_z | Filling | q | E | Full-state S | Loop purity error | Energy variance |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| armchair | 20 | 2 | -8.0 | 0.100000000 | 0 | -2.73205080756888 | 0.604736558362 | 6.44e-16 | 3.55e-15 |
| armchair | 20 | 2 | -8.0 | 0.100000000 | 1 | -2.82842712474619 | 1.124670289238 | 3.33e-16 | 3.55e-15 |
| zigzag | 19 | 2 | -7.5 | 0.105263158 | 0 | -2.73463474709466 | 0.506625559648 | 1.8e-17 | 1.78e-15 |
| zigzag | 19 | 2 | -7.5 | 0.105263158 | 1 | -2.82842712474619 | -0.000000000000 | 2.48e-16 | 0 |

All number sectors were fixed by U1 MPS QNs; winding was isolated by the exact projector penalty, with the final microscopic loop measured explicitly. Energies agree with original-basis sector ED to less than 1e-14. Four states passed the unchanged driver convergence criteria. Normal CI does not run these optimizations.

Degeneracies remain important: independent rotated-block ED finds two ground states in each zigzag N2 winding block, four in armchair q0, and one in armchair q1/q2. Entropies of degenerate candidates are not universal values and have not been identified as MES. No gamma or deconfinement claim follows.

Machine-readable values, errors, maximum dimensions, runtime, source/result paths, seeds and launch commit hashes are preserved in `sector_penalty_comparison.json`, `sector_penalty_states.csv`, and the individual result JSON files. Immutable physical-spin MPS/checkpoint payloads remain local and excluded from Git.

Library construction checks passed 36 penalty checks, 24 translated circumference checks, three local-loop checks and three translated legacy-entry checks. The direct QN basis implementation is the next active validation. Existing larger optimizations/checkpoints are preserved.
