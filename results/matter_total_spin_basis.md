# Exact total-spin decomposition and the production winding QN

No state or spin sector is truncated. Define |b> as the computational ket
with only matter spin b+1 up, and |bar b> with only spin b+1 down, b=0,1,2.
Choose C|b>=|b+1 mod 3>, and likewise for the hole kets. The coupled basis is

\[
|3/2,-3/2\rangle=|\downarrow\downarrow\downarrow\rangle,\quad
|3/2,-1/2\rangle=\frac{1}{\sqrt3}\sum_b|b\rangle,\quad
|3/2,+1/2\rangle=\frac{1}{\sqrt3}\sum_b|\bar b\rangle,\quad
|3/2,+3/2\rangle=|\uparrow\uparrow\uparrow\rangle.
\]

The multiplicity k=1,2 labels the two doublets:

\[
|1/2,-1/2;k\rangle=\frac{1}{\sqrt3}\sum_b\omega^{-kb}|b\rangle,\qquad
|1/2,+1/2;k\rangle=-\frac{1}{\sqrt3}\sum_b\omega^{-kb}|\bar b\rangle.
\]

Their cycle eigenvalues are omega^k at both magnetizations. The minus sign
fixes the conventional positive total-spin raising matrix element. In column
order (quartet at increasing m, k=1 doublet at increasing m, k=2 doublet),
and computational row order
(ddd,udd,dud,uud,ddu,udu,duu,uuu), the exact unitary is

\[
U=\frac{1}{\sqrt3}\begin{pmatrix}
\sqrt3&0&0&0&0&0&0&0\\
0&1&0&0&1&0&1&0\\
0&1&0&0&\bar\omega&0&\omega&0\\
0&0&1&0&0&-\omega&0&-\bar\omega\\
0&1&0&0&\omega&0&\bar\omega&0\\
0&0&1&0&0&-\bar\omega&0&-\omega\\
0&0&1&0&0&-1&0&-1\\
0&0&0&\sqrt3&0&0&0&0
\end{pmatrix}.
\]

Computational coordinates transform with U dagger. Consequently
U† C U = I4 direct-sum omega I2 direct-sum conjugate(omega) I2.
The complete matrix, labels and conventions are machine-readable in
`geometry/matter_total_spin_basis.json`.

For the production encoder B, U=B[:,(0,1,6,7,2,3,4,5)] times
diag(1,1,1,1,1,-1,1,-1), with zero-based column indices. Thus the existing
individual two-state mode sites already encode exactly these spin and
multiplicity states. N equals the number of occupied mode bits and
k=sum_a a*n_a mod 3. The microscopic winding generator has charge
sum_v(c_v*N_v+s_v*k_v)+sum_e p_e*n_e mod 3, where c_v=p(leg1) and
s_v=p(leg2)-p(leg1). The multiplicity label alone is insufficient: the
uniform matter number phases and physical shared gauge spins are required.
The existing onsite QN weights c_v+s_v*a implement this complete charge.

Local S is not conserved by the microscopic Hamiltonian. Matrix elements
mix the quartet and doublets; using S as a block constraint or dropping
either doublet would change the model. Only physical total number and the
specified exact winding charge constrain the production MPS.

Independent verification covers the entire eight-state matter space and all
64 states of the full single star, including exact total-spin and ladder
matrices, C, the transformed exchange operators, Hermiticity, and all energies.
Full Hamiltonian matrix error is 5.79e-16, spectral error 2.89e-15; the
signed permutation relates the coupled and production Hamiltonians with
error 2.43e-15. Local-operator equality extends to the same shared-edge lattice
by linearity, without duplicating gauge spins. Julia regression checks the
explicit unitary against the independent matrix and its exact QN encoding.

This is the already validated production representation, reindexed and
rephased, so its existing ED/projector/DMRG comparisons need no rerun.
`sector_penalty_comparison.json` records both methods on the smallest cylinder
fixtures, and `sector_preserving_physical_filling_audit.md` records the
converged physical-filling zigzag q0 state against independent ED. The q1/q2
physical-filling optimizations continue separately. This representation check
is classified REJECTED / DIAGNOSTIC ONLY for phase interpretation; no entropy
intercept or phase conclusion follows from the basis algebra alone.
