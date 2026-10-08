# Bounded physics extraction: exact microscopic WXY model

Analysis date: 2026-10-07. Source revision: `0cf8270ddab456b9487c360e5fcaed281baadf44`, with pre-existing scientific outputs preserved. This extraction finished within the 48-hour budget, using measurements and small exact symmetry-block calculations rather than new production optimizations. Energies use J=1; entropies use natural logarithms. Hamiltonian, numerical tolerances, sectors and production targets were not changed.

The strongest results are an analytical **conditional exclusion** and an exact finite-system explanation: a half-filled, gapped, primitive-translation-preserving, U(1)-preserving phase cannot have *ordinary* D(Z3) as its entire intrinsic topological order. The model's actual thermodynamic filling and symmetry realization are not established. All four tiny-torus ground states form one symmetry-enforced multiplet caused by the periodic repeated-edge geometry. Existing cylinders do not provide a converged topological-entropy inference.

## 1. Filling constraint and its limits

A honeycomb primitive cell contains six local matter spins and three physical shared gauge spins, hence nine spin-1/2 sites. The conserved physical integer charge is n=S^z+1/2 on each site. Half filling means ν=9/2 charges per primitive cell, with fractional part 1/2. Counting a shared edge twice would invalidate this statement; the saved geometry assigns it once.

Ordinary bosonic D(Z3) has nine Abelian anyons A=Z3×Z3, each satisfying a^3=1. Their fractional U(1) charges obey 3Q_a=0 modulo integers, so Q_a can only be 0, 1/3 or 2/3. No anyon can carry fractional charge 1/2. Connected onsite U(1) cannot permute discrete anyon types. The filling/LSM anomaly requires a fractional-charge spinon matching ν modulo one in a gapped phase with unbroken primitive translations and U(1). Thus ordinary D(Z3) cannot match the half-filling anomaly. The continuous-symmetry result survives translation permutations of anyons; this is explicitly discussed in Section VII A of [Cheng et al., PRX 6, 041068](https://arxiv.org/abs/1511.02263), with the cylinder argument developed by [Zaletel and Vishwanath](https://arxiv.org/abs/1410.2894).

For clarity, the obstruction is not an assumption that translations act trivially on anyons. Translations commute with physical U(1), hence preserve fractional charges. The anyon v associated with adiabatic 2π U(1) flux is translation invariant and has v^3=1. On a gapped, purely topological ground space, three flux insertions act as a scalar, up to exponentially small finite-size corrections. The flux-threading algebra is T_x F T_x^-1=exp(2πiνL_y)F. At half filling and odd L_y, cubing this relation gives a minus sign, contradicting scalar F^3. This uses the standard SET/adiabatic flux-threading assumptions, including the appropriate treatment of anyon-permuting cylinder closures; it is not a bare ground-degeneracy determinant argument on an arbitrary twisted torus. The equivalent cylinder charge argument says that primitive translation shifts fractional cut charge by one half on odd circumference, whereas changes between D(Z3) topological sectors can supply only thirds. See also the flux-insertion discussion in [Bultinck and Cheng](https://arxiv.org/abs/1808.00324).

An exact finite-arithmetic check in `physics48_filling.json` enumerates all four twist-preserving braided automorphisms of D(Z3), all sixteen commuting translation-permutation pairs, and compatible U(1) charge characters. None supplies half charge. This verifies the model-specific anyon arithmetic; the LSM/SET implication is analytical, not a numerical theorem proved by that enumeration.

**Precisely excluded:** a gapped phase at half filling with unbroken physical U(1), unbroken primitive translations, ordinary commuting translations (no external magnetic-translation algebra), and ordinary D(Z3) as the full topological order, without extra symmetry-breaking ground degeneracy. This does not prove that the actual WXY ground state satisfies those premises.

There is an exact *unitary* particle-hole symmetry: flip all physical spins and exchange matter species 2 and 3 at every star. Spin flip interchanges the two hopping orientations; the species exchange conjugates the Fourier W matrix back to itself. This commutes with H and sends N_up to N_spins−N_up. Its full N_up=9 torus matrix commutator is zero. An unbroken particle-hole symmetry in a pure thermodynamic state fixes mean half filling. Symmetry of the Hamiltonian alone does not pin the filling of a pure ground state: particle-hole symmetry can break, and finite fixed-number calculations do not resolve that possibility.

The existing number bound, reused without repeating its calculation, only establishes the necessary density window [0.1888047466362023, 0.8111952533637977]. It does not select half filling. For a symmetric ordinary D(Z3) SET the filling arithmetic requires 3ν integer, giving physical density ρ=k/27; the closest allowed densities to half are 13/27 and 14/27. These are necessary possibilities, not predicted ground densities.

Alternatives remain: different filling (and broken particle-hole symmetry for a pure state), primitive translation breaking (even cell enlargement removes this particular anomaly), U(1) breaking, gaplessness, or other/additional topological order capable of matching the half-charge anomaly. Ordinary D(Z3) coexisting with translation breaking is not excluded. The previous `primitive_filling_audit.json` and `number_bound_audit.json` remain intact.

## 2. Exact classification of the four torus ground states

The saved four eigenvectors in `torus_groundspace.npz` were reused. No ground-state Lanczos calculation was repeated. Their Gram error is 2.29×10^-14 and maximum H-eigenvector residual is 1.69×10^-14. The N_up=9 basis has dimension 48,620 and contains twelve matter plus six physical gauge spins.

Zero-based local legs for vertices A0,B0,A1,B1 are

```
[[0,1,2], [0,4,2], [3,4,5], [3,1,5]]
```

This is the periodic honeycomb cell, including parallel connections, not K4. For gauge exponents p_e modulo three with sum zero at each vertex, define c_v=p_(leg0), r_v=p_(leg1)−c_v. The exact microscopic CGS action permutes matter a→a+r_v and multiplies a computational state by ω^(Σ_v c_v N_matter,v + Σ_e p_e n_e). All matter action is retained.

| Operator | Gauge exponent vector p |
|---|---|
| Y0, short circumference winding at pair 0 | (1,0,2,0,0,0) |
| Y1, short circumference winding at pair 1 | (0,0,0,1,0,2) |
| X, longitudinal winding | (1,2,0,1,2,0) |
| P_hex, folded contractible hexagon | (1,0,2,2,0,1) |

In this number sector Y0 Y1 X=I and P_hex=Y0 Y1†. The first relation uses total charge N_up=9; it is not a sector-independent operator identity. The cycle operators commute. Y0,Y1 supply a maximal commuting set on the four-dimensional ground space. Charges q mean eigenvalue ω^q.

| Ground state | q(Y0) | q(Y1) | q(X) | q(P_hex) | Translation T_x |
|---|---:|---:|---:|---:|---|
| 0 | 2 | 2 | 2 | 0 | +1 |
| 1 | 2 | 1 | 0 | 1 | exchanges with state 2 |
| 2 | 1 | 2 | 0 | 2 | exchanges with state 1 |
| 3 | 1 | 1 | 1 | 0 | +1 |

All four energies are −7.364534839426. T_y is the identity on this one-cell circumference quotient. Particle-hole F interchanges both charges; the individual charge basis states are not F eigenstates. F's zero expectation in those states must not be described as an eigenvalue.

The crucial additional exact symmetries are R0 and R1. R_k exchanges the two parallel gauge edges at pair k. At both adjacent stars it maps matter a→−a modulo three and multiplies by ω^(2Σ_a a n_a). The phase is required by the microscopic W matrix. These satisfy R_k²=I, R0R1=R1R0, R0Y0R0=Y0† and R1Y1R1=Y1†, with each reflection commuting with the other cycle. Thus the exact symmetry contains D3×D3 (each D3 is the order-six dihedral group). R0 flips the first nonzero charge and R1 flips the second, connecting all four states. The ground space is one irreducible two-dimensional doublet tensor another doublet: **all fourfold degeneracy is explained**, with no additional degeneracy left over. A different commuting basis can diagonalize R0,R1 with all four ± parity combinations, but cannot also assign definite Y0,Y1 charges.

These reflections exploit parallel edges with identical endpoints and do not extend as these local operators to a generic wider honeycomb lattice. Consequently this degeneracy is not evidence for thermodynamic topological protection. The ground states also have different contractible P_hex charges. That is incompatible with interpreting these four particular states as an already demonstrated locally indistinguishable topological ground manifold. It does not exclude a different thermodynamic phase.

All proposed operators were checked against the complete microscopic Hamiltonian in the full N_up=9 block, not merely within the ground space. Maximum commutator entry is 3.56×10^-16, ground-space invariance residual is at most 6.80×10^-13, and exact group-relation errors are at most 1.1×10^-15. Full 4×4 complex operator matrices, restricted commutators and symmetry connection matrices are in `physics48_torus.json`.

### Sector excitations

Saved excited eigenvectors were unavailable. Only the missing excited levels were therefore computed in small exact CGS orbit blocks, with cached ground vectors deflated rather than recomputed. Lanczos tolerance was 10^-12, seed 7310 with charge-dependent offsets. No tolerance was loosened. Residuals are at most 7.36×10^-13.

| Charge orbit (Y0,Y1) | Dimension per sector | Lowest energy | Next distinct energy |
|---|---:|---:|---:|
| (0,0) | 5636 | −7.342883075944 | −7.337288472473 |
| (±1,0), (0,±1) | 5454 | −7.296046028855 | −7.289970027036 |
| (±1,±1) | 5292 | −7.364534839426 | −7.136881951623 |

Dimensions sum to 48620. Symmetry relates the four members of each nontrivial orbit. This reproduces multiplicities 4,1,1,4 in the corrected first ten levels, and the next four at −7.289970027036. The global finite-cell gap is 0.021651763482; the first excitation within each ground cycle sector is 0.227652887803 above it. Neither is a thermodynamic gap or a deconfinement measurement. `results/torus_hard_gate_audit.md` is preserved as the historical validation record.

### Tiny-torus spatial entropy

For the prior full-state cut, left physical bits are {0,1,2,3,4,5,12,14,16}. Y0,R0 act entirely on the left; Y1,R1 entirely on the right. Each definite cycle-charge state has S_base=1.275330241624596. Arbitrary ground-space coefficients form a two-by-two matrix c in the reflection-consistent orbit basis. Direct Schmidt decomposition gives

```
S = S_base + entropy(singular_values(c)^2),
S_base <= S <= S_base + ln(2).
```

A logical Bell superposition has S=1.968477422184540. Product, Bell and four seeded random coefficient matrices verify the formula to 1.46×10^-13. This explains the entropy variation as ordinary entanglement between the two encoded symmetry doublets. It is not thermodynamic TEE, and the earlier entropy outputs are retained for provenance.

## 3. Existing cylinder states: measurements and rejection reasons

Twelve width-two/three checkpoints were measured without optimization: L=4 at χ256 and χ512, and L=6 at χ256, for both families. Checkpoint hashes were checked before and after measurement. All finite fillings below are fixed by saved total-U(1) MPS sectors, not chosen by minimizing over sectors. Open boundaries have extra dangling gauge spins; physical spin counts therefore differ from 9wL.

Geometry uses the existing explicit `geometry/{family}_L{L}_w{w}_star.json` and CGS-cycle tables. Zigzag wraps w a2, with longitudinal row t=x; armchair wraps w(a1+a2), with t=x−y. Both have open longitudinal boundaries and one spin per physical shared edge. In nearest-neighbor-distance units, circumferences are sqrt(3)w and 3w, respectively. Individual physical sites use the saved star ordering; the full-state cut is the saved spatial partition, including matter and a fixed shared-edge assignment. Different families are not pooled into one fit.

### Best saved L=4 states (χ512)

| Family, width | Spins | N_up; S^z; filling | E | Full S | Dominant winding q; weight | 1−|〈U〉|² | Central projected residual | Last-pass max truncation |
|---|---:|---|---:|---:|---|---:|---:|---:|
| Zigzag 2 | 74 | 37; 0; 0.5 | −29.9351369562 | 2.118235578 | 0; 0.999422963 | 0.00173036 | 0.00378473 | 0.000306285 |
| Armchair 2 | 76 | 38; 0; 0.5 | −30.6349020760 | 3.084902398 | 2; 0.999702389 | 0.000892634 | 0.00556773 | 0.000159572 |
| Zigzag 3 | 111 | 56; +0.5; 0.504504505 | −44.7733929283 | 2.138517948 | 1; 0.999751019 | 0.000746800 | 0.00162213 | 0.0000964741 |
| Armchair 3 | 114 | 57; 0; 0.5 | −45.6643879396 | 3.420974201 | 0; 0.999509180 | 0.00147192 | 0.00558853 | 0.000150470 |

The central winding MPO includes the exact matter permutations and phases. Its conjugate/U² identity agrees to 10^-10. These states overwhelmingly favor one exact CGS charge, but they do not satisfy the previous 10^-10 purity criterion. That is a small quantified admixture, not proof of a large sector cat. A pure microscopic CGS charge would itself still require a separate identification with the topological flux/MES relevant to TEE.

Matched continuation branches retain their dominant winding charge between χ256 and χ512, but show the following changes:

| Family, width | ΔS, χ256→512 | ΔE/vertex, χ256→512 | Last same-χ refinement ΔS |
|---|---:|---:|---:|
| Zigzag 2 | +0.490242 | −0.00192690 | 0.00559608 |
| Armchair 2 | +0.332020 | −0.00629446 | 0.124615 |
| Zigzag 3 | +0.249425 | −0.00268151 | 0.0292151 |
| Armchair 3 | +0.443820 | −0.0205728 | 0.000678348 |

All exceed the existing bond-entropy target 0.001, refinement-entropy target 0.0001 and bond-energy-per-vertex target 10^-6. Truncation errors exceed the existing 10^-8 target by four orders of magnitude. The central two-site residual is a retained-variational-space stationarity diagnostic, not a full H-eigenvector residual or a ground-state certificate. No new acceptance threshold was invented for it.

### Length information, χ256

| Family, width, L | N_up; S^z; filling | Full S | Dominant q | Winding variance |
|---|---|---:|---:|---:|
| Zigzag 2, 6 | 55; 0; 0.5 | 1.592340488 | 1 | 0.00103123 |
| Armchair 2, 6 | 56; 0; 0.5 | 2.439571423 | 2 | 0.000763108 |
| Zigzag 3, 6 | 83; +0.5; 0.503030303 | 1.874589667 | 0 | 0.00219715 |
| Armchair 3, 6 | 84; 0; 0.5 | 3.018535364 | 1 | 0.00246302 |

Most L4/L6 comparisons change dominant charge and cannot be used as fixed-sector length convergence. Armchair width-two has the same nominal half filling and dominant q=2, but its χ256 entropy changes by −0.313311. Initializations are independent and neither sector is pure, so this is only suggestive length sensitivity, not a controlled convergence measurement.

For physical density profiles, each row includes its six-w matter spins and three-w gauge spins allocated once to their unwrapped A owner; outside-owner dangling edges contribute to total filling but not these row averages. The central L4 row densities (t=1,2) are: zigzag w2 (0.497852,0.482381), armchair w2 (0.516488,0.484091), zigzag w3 (0.498249,0.517897), armchair w3 (0.487647,0.513190). The centers are not demonstrably homogeneous. These short, unconverged, boundary-affected profiles cannot distinguish spontaneous translation breaking from boundaries or variational error. Full physical-site profiles and length-six row averages are saved in supporting JSON.

Existing finite number scans do not establish the bulk filling: their energies are variational upper bounds, convergence errors compete with differences between number sectors, and comparing distinct χ/initialization branches is not a reliable sector ordering. No scan was restarted.

### Infinite saved states and transfer interpretation

The saved armchair width-two nominal χ64 state actually has maximum bond dimension nine and a one-dimensional cell-wrap transfer space. It has enforced background N_up=18, S^z=0 per 36-spin/two-slice cell, filling 1/2. Its source solver residual is 3.40×10^-9 and source canonical error is 3.40×10^-9. A measurement-only recanonicalization of the same AL state was justified by a complete dense transfer-space certificate: the sole eigenvalue is 0.9999999999999997+7.18×10^-18 i, residual zero at numerical precision. No optimization was performed.

Exact remeasurement finds pure winding charges alternating q=2 and q=1 on consecutive slices, variances below 5.0×10^-15, and contractible hexagon q=1 for the tested representatives. The row fillings are 8/18=4/9 and 10/18=5/9. This saved state explicitly breaks the one-slice primitive translation while preserving half filling on average. It is **not** a candidate satisfying the translation premise of the half-filling no-go, and its purity does not certify it as the global ground state. Its energy/vertex is only −1.436915904, distinctly different from the other saved branches. Its essentially zero slice-cut entropy and scalar wrap transfer describe a cell-product variational state. There is no resolved subleading eigenvalue: a recorded ξ=0 convention is not evidence for a microscopic spectral gap. Within-cell correlations can still be nonzero.

The zigzag width-two χ64 sequential and parallel branches record, respectively, E/vertex −1.816335670 and −1.823540911, slice entropy averages 0.676018 and 0.765264, and ξ_cell≈0.149258 and 0.162274. Their solver residuals are 0.00428738 and 0.00319015; canonical errors are 0.00428732 and 0.00359080. Both zigzag states use an unrestricted-number ansatz: N_up and S^z are expectations, not enforced cell quantum numbers. Sequential: N_up≈17.999998618, S^z≈−0.000001382 per 36-spin cell, filling≈0.499999962. Parallel: N_up≈18.000095697, S^z≈0.000095697, filling≈0.500002658. Neither establishes variational selection of the global ground filling. These are unconverged variational transfer spectra, not trustworthy physical bulk correlation lengths or a bond-converged sequence. Leading real/imaginary eigenvalues and eigenpair residuals are preserved in `physics48_summary.json`; even a small transfer-eigenpair residual cannot repair a nonstationary or improperly canonicalized variational state. A transfer length alone never certifies the microscopic Hamiltonian gap.

### Reliability and admissibility

- **Reliable as saved-state properties:** finite fixed N_up/S^z, full-state cut entropies, density profiles, untruncated loop expectations, quantified sector admixtures, checkpoint identity; exact torus matrices and spectrum; the certified scalar transfer spectrum of the recanonicalized armchair cell-product state.
- **Suggestive but unproven physically:** charge-dominated finite states, density modulation/boundary sensitivity, the period-two infinite candidate, finite-sector energy preferences.
- **Unconverged or rejected for phase inference:** all twelve finite cylinder candidates, both zigzag infinite transfer-length branches, thermodynamic filling selection, microscopic bulk gap, intrinsic topological flux/MES identification, sector splitting extrapolation and any gamma intercept.

There are **zero accepted fit points**. No gamma fit, uncertainty, or comparison to ln(3)=1.0986122886681098 is scientifically admissible. Raw entropies are retained separately from this decision. Two narrow circumferences plus the recorded drifts cannot establish a 2D phase. No independently validated microscopic string/deconfinement calculation exists in this extraction, and none was invented.

## 4. Scientific decision

1. **Established rigorously:** the conditional half-filled symmetric ordinary-D(Z3) exclusion, exact particle-hole symmetry without a ground-filling theorem, and the finite torus D3×D3 explanation. Numerically, the saved torus vectors, group actions and corrected sector spectrum are validated at their stated residuals.
2. **Strongly suggested but unproven:** the fourfold tiny-torus degeneracy is a geometry/symmetry effect rather than independent evidence for topology. Existing optimizations explore materially different, sometimes translation-modulated branches; their phase cannot be inferred from one branch or an apparent intercept. No particular thermodynamic symmetry-breaking pattern is established.
3. **Evidence concerning intrinsic Z3 order:** there is no reliable positive numerical evidence. There is a firm analytical obstruction to *ordinary symmetric D(Z3) at half filling*, conditional on actual filling and symmetry assumptions. There is no unconditional exclusion of intrinsic D(Z3) coexisting with broken symmetry or occurring at another filling, nor of more complicated topological order.
4. **Single most important unknown:** which filling and primitive-translation realization belongs to the true zero-field thermodynamic ground state. This determines whether the analytical exclusion applies; the present fixed-sector/period-restricted variational branches do not answer it.
5. **Can the current approach realistically decide?** ITensor DMRG/VUMPS can in principle address it, but the observed Mac-scale data do not support a credible near-term TEE extrapolation: entropy drift remains 0.25–0.49 upon χ doubling, purity fails, and several infinite branches are nonstationary. Merely extending the current queue is not a justified promise of a definitive 2D answer. A clean period-two stationary candidate is not a global minimization proof.
6. **Cheapest decisive next step:** pursue a model-specific analytical theorem or sharp filling-dependent energy bound fixing the actual ground filling *and* translation realization. If half filling and unbroken translations/U(1) can be established, ordinary D(Z3) is excluded without larger cylinders. If primitive translation breaking is established, this no-go no longer decides the intrinsic order. No currently identified inexpensive finite calculation guarantees an unconditional answer; small-cluster ED or more unconverged entropy points should not be advertised as decisive. A matched-sector, matched-period variational comparison would be useful only after demonstrating credible stationarity/convergence, which the existing wider-cylinder data presently lack.

## Preservation, validation and reproduction

The running production process PID 74718 was paused with SIGSTOP after its latest checkpoint hash and actual resume-loader validation succeeded. Its state and parent were retained. The lightweight watcher was given an intentional user-authorized pause marker so it cannot launch the next expensive job. The verified χ16 iteration-30 checkpoint and unfinished χ128 target are recorded in `physics48_production_pause.json`; uncheckpointed in-memory progress is preserved. Production remains paused at this report's completion because this request prohibits new large-scale optimization. No queue redesign was made.

Julia 1.10.10 and the existing Project.toml/Manifest.toml were used with one BLAS thread; analysis used the existing Python 3.14 environment and requirements. Measurement regressions passed: 24 exact/translated CGS MPO checks, 36 spatial-partition checks and five stationarity-fixture checks. The torus symmetry-block levels reproduce the corrected benchmark; new entropy decompositions have independent exact coefficient-matrix checks. The legacy L6 checkpoint aliases initially lacked a manifest pointer; the loader was corrected to accept the existing valid alias, and measurements succeeded without optimization. Error records remain as provenance.

Run from repository root with the preserved environment:

```
python scripts/physics48_filling.py
python scripts/physics48_torus.py
python scripts/physics48_torus_excitations.py
python scripts/physics48_torus_entanglement.py
OPENBLAS_NUM_THREADS=1 julia --project=. scripts/physics48_cylinders.jl
OPENBLAS_NUM_THREADS=1 julia --project=. scripts/reaudit_infinite_loops.jl results/infinite_armchair_w2_chi64_qn_star_active.jls
python scripts/physics48_summary.py
```

The measurement commands require the local saved binary checkpoints, deliberately excluded from ordinary Git. The new small JSON results preserve all numerical conclusions, source hashes, settings and data needed to inspect this audit without those binaries. Original simulation scripts and environment remain the route to regenerating checkpoints from scratch. The GitHub report must not be construed as distributing the original saved wavefunctions. New torus joint/excited vectors are also local ignored binaries. No essential new inference exists only in terminal output.
