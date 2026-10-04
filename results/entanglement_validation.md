# Entanglement and transfer validation

The same `entropy_at` routine used for WXY finite cylinders passes independent
exact Schmidt-spectrum tests (`test/entanglement.jl`): a product state has
S=0; a normalized two-term state with weights 0.3 and 0.7 has those Schmidt
probabilities and S=-0.3 ln(0.3)-0.7 ln(0.7); equal weights give S=ln(2).
All four assertions pass at 1e-14 absolute tolerance. Independently, the
19-spin zigzag L=2,w=1 full spatial bipartition gives
S_ED=1.3553089306439614 and S_ITensor=1.3553089188696528. The 20-spin armchair
groundspace is degenerate, so the entropy of one arbitrary ED eigenvector is
not a unique reference for the selected DMRG state.

`test/transfer.jl` constructs an analytic two-dimensional-bond channel from
A_up=diag(sqrt(0.8),sqrt(0.2)) and
A_down=[[0,sqrt(0.8)],[sqrt(0.2),0]]. Its eigenvalues are 1,0.8,0,0, so
xi=-1/ln(0.8)=4.48142011772455 in one-site cells. Library TransferMatrix and
KrylovKit recover the leading eigenvalues, xi, and residuals within the
recorded 1e-10/1e-9 tolerances.

A separate GHZ transfer test has eigenvalues 1,1,0,0. A scalar Krylov space
can hide the repeated fixed point, just as scalar ARPACK hid torus energy
multiplicities. `transfer_spectrum` checks independent dominant fixed-point
starts and their Gram matrix; the GHZ test detects rank two and infinite xi.
All seven transfer assertions pass. Near-unit ratios unresolved relative to
residuals are explicitly flagged; a tiny computed nonzero subleading value
below numerical resolution is not a precise physical correlation length.

The audited checkpointing VUMPS driver calls the package's existing VUMPS
iteration, rather than implementing tensor-network operations. Its energy
and center normalization agree with the standard package driver in the
one-iteration regression. The finite checkpoint observer preserves the
single-star benchmark and resumes the hashed saved state without changing
energy or entropy. Large state payloads remain local; small manifests and
validation logs are preserved.

These tests validate the measurement pipeline, not WXY topological order.
That claim still requires ground-state, bond-dimension, length, circumference,
and sector convergence for both cylinder families.
