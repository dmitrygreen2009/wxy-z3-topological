"""Decide the torus hard gate using full spectra and independent certificates."""
import json,gzip,csv,pathlib
import numpy as np

independent=json.load(open('results/torus_independent.json'))
audit=json.load(open('results/torus_convention_audit.json'))
solvers=json.load(open('results/torus_solver_audit.json'))
expected=np.array([-7.364534839425767]*4+[-7.342883075944,-7.337288472473]+[-7.296046028855]*4+[-7.289970027036]*4)
assert max(abs(np.array(independent['energies'][:14])-expected))<1e-9
assert independent['converged']>=16
assert independent['gram_error']<1e-8 and max(independent['residuals'])<1e-9
assert audit['enumerated_patterns']==62208 and audit['matching_patterns']==0
assert all(g['ground_multiplicity']==4 and g['all_basis_popcounts_are_9'] for g in audit['translation_quotients'])
assert all(max(abs(np.array(g['energies'][:14])-expected))<1e-9 for g in audit['translation_quotients'])
with gzip.open('results/torus_all_endpoint_conventions.csv.gz','rt') as f:
    rows=list(csv.DictReader(f))
assert len(rows)==62208
assert all(max(abs(np.array([float(r[f'E{k}']) for k in range(8)])-expected[:8]))<1e-9 for r in rows)
assert sum(p['dimension'] for p in solvers['parity_blocks'])==48620
assert all(sum(abs(np.array(p['energies'])-expected[0])<1e-9)==2 for p in solvers['parity_blocks'])
result={'gate_condition_satisfied':2,'supplied_full_first8_reproduced_as_complete_spectrum':False,
    'supplied_full_first8_independently_refuted':True,'corrected_first8':expected[:8].tolist(),
    'ground_multiplicity':4,'independent_solver':'Julia eight-vector BlockLanczos',
    'endpoint_cases_checked':len(rows),'scalar_arpack_can_return_supplied_incomplete_list':any(r['matches_supplied_first8'] for r in solvers['scalar_arpack']),
    'larger_system_scientific_conclusions':'Not established; original scaling study resumed after user acceptance of corrected benchmark.'}
pathlib.Path('results/torus_gate.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
