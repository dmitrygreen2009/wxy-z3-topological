"""Summarize exact microscopic loops; archived triplets are excluded."""
import hashlib, json, pathlib, subprocess
rows=[]
for p in sorted(pathlib.Path('results').glob('infinite_*_circumference_cycles.json')):
    d=json.loads(p.read_text())
    if 'noncontractible_cycle_measurements' not in d:continue
    rows.append({'source':str(p),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
        'family':d['family'],'width':d['width'],
        'source_solver_residual':d.get('source_solver_residual'),
        'source_canonical_error':d.get('source_canonical_error'),
        'circumference_charge_probabilities':[r['charge_probabilities'] for r in d['noncontractible_cycle_measurements']],
        'ground_state_or_2D_order_certified':False,
        'translated_operator_reaudit_status':'Pending saved-state remeasurement; no flux/MES certificate'})
out={'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'records':rows,
     'excluded_observables':['matter-triplet one-point functions','matter-triplet correlations'],
     'interpretation':'Exact microscopic loop measurements of saved variational states. Loop purity does not certify a topological flux sector or MES. Translated-indexing audits and sector-preserving optimization remain required.'}
pathlib.Path('results/infinite_microscopic_observables_summary.json').write_text(json.dumps(out,indent=2)+'\n')
