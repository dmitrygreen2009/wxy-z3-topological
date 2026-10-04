# Global number sectors on microscopic open cylinders

The benchmarks specify their number blocks explicitly. Their successful
validation does not imply that half filling minimizes every open-cylinder
Hamiltonian. The complete independent Python ED scan in
`small_cylinder_number_sectors.json` checks every N_up for both L2 width-one
geometries, including residuals and particle-hole partner energies.

| Family | Physical spins | Minimizing N_up | Global E0 | Half-filled E0 |
|---|---:|---|---:|---:|
| Zigzag | 19 | 9, 10 | -7.696107978815068 | -7.696107978815051 |
| Armchair | 20 | 8, 12 | -8.037741645057586 | -7.994302046793339 |

Maximum scan residuals are 6.19e-11 and 7.94e-11. Scalar k=1 certifies
energies by number sector, not complete degeneracies. The exact energy
advantage of the armchair off-half-filled sectors is approximately 0.04343960.
Independent Julia eight-vector block ED reproduces both minimizing energies
with one ground state in each sector and residuals below 1e-12. Together with
the complete number scan this certifies a two-state global groundspace. The
within-sector excitation gap is 0.5424605484 in each minimizing sector.
Microscopic ITensor DMRG at maximum bond dimension 256 (achieved 230)
reproduces the N_up=8 energy within 6.3e-11. Its full-state spatial entropy
is 1.14414415305 and direct energy variance is 5.45e-10. Quantum-number
labels are checked against the actual MPS flux. All shared edges remain one
physical spin, and the local matter spins remain in the spatial entropy.

Existing armchair half-filled cylinder points must be interpreted as
fixed-number-sector results. This tiny-cylinder observation does not imply
that the preferred filling shift remains extensive in the thermodynamic
limit: boundary effects and bulk density selection require longer cylinders.
No fit of the existing armchair entropies alone determines the global-ground-
state topological constant. The same global-sector checks must accompany
larger zigzag calculations even though this small zigzag example minimizes
at the two closest half-filled sectors.

The finite driver now accepts an explicit `nup` keyword, labels the sector in
result/checkpoint names, checks the MPS quantum-number flux, and preserves
it on continuation/resume. Charged searches include offsets of two as well
as one and three; the former cannot be omitted on this geometry. Independent
higher-dimensional sector comparisons remain variational until convergence
and competing sectors have been examined.

Reproduce the independent scan:
```sh
OPENBLAS_NUM_THREADS=1 python scripts/benchmark_number_sectors.py --cylinders
julia --project=. scripts/cylinder_block_ed.jl armchair:8 armchair:12
julia --project=. -e 'include("scripts/cylinders.jl"); run_cylinder("armchair",2,1,128;ordering="star",nup=8)'
```
