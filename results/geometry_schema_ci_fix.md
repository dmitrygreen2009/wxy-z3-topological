# Cylinder geometry schema regression (2026-10-08 UTC)

The CI failure in `scripts/strategy_geometry_audit.py` was an omitted field,
not a renamed field. The committed `geometry/armchair_L6_w3_star.json` had
legacy conventions without `physical_axial_matter_offsets_A_B` or
`physical_axial_gauge_offsets_from_owner_A_by_leg`. The current
`src/geometry.jl::export_cylinder_geometry` already emits both fields.

Restore the exporter metadata using `scripts/normalize_geometry_metadata.py`.
The physical graph, endpoint-leg table, matter spins, individual-site MPS
ordering, and spatial cut are unchanged. Armchair offsets are `[0,0]` for
matter and `[0,-1/2,1/2]` for gauge legs; zigzag offsets remain `[0,1/3]`
and `[1/6,-1/3,1/6]`. No Hamiltonian or numerical tolerance changes.

The auditor still requires these fields and independently verifies physical
cut membership using integer sixth-slice coordinates, physical shared-edge
counts, periodic endpoint translations, and closed microscopic winding paths.
Missing metadata now gives an explicit assertion with the offending path and
export command instead of an unexplained `KeyError`.

Permanent regressions check all stored manifests and freshly generated
zigzag/armchair exports. Negative fixtures must reject a missing offset,
an incorrect offset, and an incorrect spatial cut.

Local affected checks passed: 41 finite-cylinder records, 40 infinite-cell
records, 221 winding cycles, four fresh-export assertions, and all three
negative fixtures. Validation ran in an isolated temporary checkout to avoid
overwriting scientific results or checkpoints. The first 29 automatic push-CI
steps passed. All non-optimization tests in the final Julia suite also passed. The final
Julia suite contains existing DMRG/VUMPS optimization
fixtures; permission to run those is pending under the user's explicit
no-optimization instruction. Production remains paused. Full CI and push
must not be reported as completed until this conflict is resolved.
