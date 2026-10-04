"""Prevent archival matter-triplet diagnostics from reentering core drivers."""
from pathlib import Path
import json
for name in ('observables.jl','infinite_observables.jl','infinite_circumference_cycles.jl'):
    text=(Path('scripts')/name).read_text()
    assert 'infinite_triplet_' not in text,name
    assert 'triplet_correlation(' not in text,name
    assert 'matter_charge_three_correlations' not in text,name
r=json.loads(Path('results/scientific_strategy_reaudit.json').read_text())
assert 'one_unrestricted_variational_ansatz_breaks_physical_U1' not in json.dumps(r)
for f in ('scientific_strategy_reaudit.md','scientific_strategy_reaudit.json','circumference_operator_validation.md'):
    assert (Path('results/archive/triplet_provenance')/('before_core_removal_'+f)).exists()
assert Path('results/torus_hard_gate_audit.md').exists()
print('Core diagnostic exclusions and archival provenance passed')
