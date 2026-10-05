# Exact real representation of the microscopic exchange

This representation follows the validated eight-state matter basis B in
`geometry/matter_charge_basis.json`. No local state is removed. In the rotated
computational basis define the three individual spin modes i,j,k, where i is
the original local leg label and j,k are the other two modes. Then

    A_i = B† (sum_a W_ai μ_a^-) B
        = a S_i^- + b S_i^- Sz_j Sz_k + d S_i^+ S_j^- S_k^-
    a = (1 - 1/sqrt(3))/2
    b = 2(1 + 1/sqrt(3))
    d = 2/sqrt(3).

The coefficient of S_i^- is one when j,k have equal occupations and
-1/sqrt(3) when their occupations differ. The last term connects the state
with mode i down and j,k up to mode i up and j,k down with amplitude
2/sqrt(3). These account for all five nonzero matrix units of A_i. The identity
was independently checked on the complete eight-dimensional matter space for
all three legs to tolerance 1e-13; errors are recorded in
`real_charge_hamiltonian_audit.json`.

The full Hamiltonian remains exactly

    H' = -sum_v,i (A_vi σ_e_i^+ + A_vi† σ_e_i^-).

There are 18 real product terms per star including Hermitian conjugates.
Physical shared edges are still represented once. This may reduce production
memory relative to a complex matrix-unit representation, and the library MPO implementation in `src/real_charge_hamiltonian.jl` has now passed its separate equivalence test in GitHub CI at commit d3e7db21. The run URL and validation scope are saved in `real_charge_hamiltonian_ci_validation.json`. Running jobs
retain their original loaded Hamiltonians. A complex checkpoint must not be
converted to a real state by taking the real part of each MPS tensor.

The onsite QN selection rule is explicit. Using zero-based mode indices, any
allowed microscopic CGS generator has gauge leg weight p_i=c+q i (mod 3) and
matter mode weight w_i=c+q i (mod 3). The first two terms have matter winding
change -w_i, canceled by +p_i on the shared gauge spin. The last term has
change w_i-w_j-w_k; adding p_i gives 3q i=0 mod 3 because i+j+k=3. All terms
lower matter number by one and raise gauge number by one. Thus each product
term preserves physical N_up and every microscopic CGS cycle charge. This
checks the exact QN structure independently of an approximate low-energy
description.

The three-mode exchange is required by the exact basis change of H. It is not
the archived S+S+S+ observable, an added interaction, or evidence for a phase.
Real coefficients also do not establish absence of sign frustration: some
exchange amplitudes change sign with other-mode occupations. No phase
conclusion is drawn from this representation.
