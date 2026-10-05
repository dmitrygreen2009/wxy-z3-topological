# Exact limitation of a two-site U1 x Z3 random initializer

The pinned ITensorMPS `random_mps` implementation randomizes a product state
with adjacent two-site unitaries conserving the supplied QNs. For two spin-half
sites with winding weights w_i and w_j, the occupation states have charge
pairs (N,W): (0,0), (1,w_i), (1,w_j), (2,w_i+w_j mod 3).
If w_i differs from w_j modulo three, every pair is distinct. Such a
two-site conserving unitary is therefore diagonal and cannot exchange their
occupations. For equal weights it can exchange one particle inside the run.
Consequently this initializer retains the total occupation of each contiguous
run of equal winding weight. These are extra initializer restrictions, not
additional microscopic conservation laws.

This follows from the exact QN block structure, not from an ad hoc physical
observable. `test/projected_number_initializer.jl` checks the diagonal form
of the official random gates for all six unequal-weight pairs. This fact
alone does not prove that later optimization noise cannot escape a particular
initialization. Energy variance and independent sector ED remain the rejection
criteria for a stalled optimizer.

The alternative in `src/projected_number_initializer.jl` first randomizes a
physical spin MPS with only the conserved physical number. It then uses the
already validated exact P_q B† checkpoint bridge, which retains all eight
matter states, to enter the desired number/winding block. The three-site local
matter basis rotation and exact projector can populate states unavailable to
the initial adjacent U1 x Z3 random circuit. Every tensor operation remains in
ITensorMPS; no tensor-network library or effective model is introduced.

Fresh finite winding jobs use this strategy after its regression gate passes;
the legacy `direct_two_site_qn_random` keyword remains available to reproduce
older runs. Resumed checkpoints retain their full amplitudes and QN sectors.
The recorded seed, physical product state, number-only random bond dimension,
projection weight, converted scalar type and compression cutoff reproduce the
initializer. Existing running jobs retain their previously loaded code and
are not interrupted. The six physical-filling geometry/charge projector
fixtures are implementation checks, not ground-state or MES evidence.
