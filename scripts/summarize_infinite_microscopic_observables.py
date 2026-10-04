"""Keep measured variational observables separate from thermodynamic conclusions."""
import hashlib, json, pathlib, subprocess
rows=[]
for p in sorted(pathlib.Path('results').glob('infinite_*_circumference_cycles.json')):
    d=json.loads(p.read_text())
    if 'matter_charge_three_correlations' not in d:
        continue
    cs=d['matter_charge_three_correlations']
    rows.append({'source':str(p),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
        'family':d['family'],'width':d['width'],
        'source_solver_residual':d['source_solver_residual'],
        'source_canonical_error':d['source_canonical_error'],
        'circumference_charge_probabilities':[r['charge_probabilities'] for r in d['noncontractible_cycle_measurements']],
        'maximum_matter_triplet_one_point_magnitude':max(abs(complex(r['one_point_real'],r['one_point_imag'])) for r in cs),
        'largest_distance_cells':max(r['distance_cells'] for r in cs),
        'largest_distance_raw_correlation_magnitudes':[abs(complex(r['correlation_real'],r['correlation_imag'])) for r in cs if r['distance_cells']==max(q['distance_cells'] for q in cs)],
        'largest_distance_connected_correlation_magnitudes':[abs(complex(r['connected_real'],r['connected_imag'])) for r in cs if r['distance_cells']==max(q['distance_cells'] for q in cs)],
        'ground_state_or_2D_order_certified':False})
out={'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'records':rows,
     'interpretation':'A nonzero CGS-invariant charge-three one-point expectation diagnoses U(1) breaking of this variational ansatz. It does not establish spontaneous order in the converged cylinder or 2D ground state. Short connected transfer correlations can coexist with its disconnected ordered plateau. Circumference charge weights do not certify a MES; separate initialization branches and convergence remain necessary.'}
pathlib.Path('results/infinite_microscopic_observables_summary.json').write_text(json.dumps(out,indent=2)+'\n')
