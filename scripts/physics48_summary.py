"""Read-only consolidation of the bounded physics extraction; no optimization."""
import json, pathlib, hashlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
def read(p): return json.loads((ROOT/p).read_text())
rows=[]
for p in sorted((ROOT/'results').glob('physics48_*_L*_w*_chi*.json')):
 if p.stem.endswith('_error'): continue
 d=json.loads(p.read_text()); src=read(d['source_result']); record=src['records'][-1]
 r={k:d[k] for k in ['family','length','width','cap','physical_spins','nup','total_sz','filling_fraction','full_state_entropy','winding_charge_probabilities','winding_unitary_variance','source_checkpoint','source_sha256','spatial_cut']}
 r.update(energy=record['energy'],dominant_q=max(range(3),key=lambda q:d['winding_charge_probabilities'][q]),projected_residual=d['bond_stationarity']['projected_residual'],maximum_truncation_error=max(record.get('sweep_max_truncation_errors',[record.get('truncation_error',float('nan'))])),entropy_change_refinement=record.get('entropy_change_refinement'),measurement_file=str(p.relative_to(ROOT)),phase_inference_eligible=False)
 g=read(f"geometry/{d['family']}_L{d['length']}_w{d['width']}_star.json"); den=d['physical_site_density_profile']; groups={t:[] for t in range(d['length'])}; vt={v['id']:v for v in g['vertices']}
 for m in g['matter_spins']: groups[vt[m['vertex_id']]['t']].append(den[m['physical_spin_id']-1])
 for e in g['shared_gauge_edges']:
  x,y=e['owner_A_unwrapped_bravais']; t=x if d['family']=='zigzag' else x-y
  if t in groups: groups[t].append(den[e['physical_spin_id']-1])
 r['owner_A_row_densities']=[{'row':t,'spins':len(vals),'density':sum(vals)/len(vals)} for t,vals in groups.items()]
 r['row_density_convention']='Matter at t and gauge allocated once to unwrapped owner A at t; outside-owner dangling edges excluded from row averages only.'
 rows.append(r)
changes=[]
for r in rows:
 if r['cap']!=512: continue
 s=next(s for s in rows if (s['family'],s['length'],s['width'],s['cap'])==(r['family'],r['length'],r['width'],256))
 changes.append({'family':r['family'],'width':r['width'],'length':r['length'],'entropy_change_256_to_512':r['full_state_entropy']-s['full_state_entropy'],'energy_change_per_vertex':(r['energy']-s['energy'])/(2*r['width']*r['length']),'dominant_q_unchanged':r['dominant_q']==s['dominant_q']})
infinite=[]
for p in ['results/infinite_armchair_w2_chi64_qn_star_active.json','results/infinite_zigzag_w2_chi64_finite_seed_star.json','results/infinite_zigzag_w2_chi64_finite_seed_star_parallel_expanded.json']:
 d=read(p); fields=['family','cap','chi','cell_spins','cell_slices','number_background','energy_per_vertex','spatial_entropy_slice_mean','solver_residual','canonical_error','xi_cells','transfer_eigenvalues_real','transfer_eigenvalues_imag','transfer_residuals','bond_dimension_details']
 profile=d['site_sz_profile_infinite_order']; sz=sum(profile); density=.5+sz/d['cell_spins']
 infinite.append({'unit_cell_sz_expectation':sz,'unit_cell_nup_expectation':sz+d['cell_spins']/2,'filling_fraction_measured':density,'filling_fixed':d['number_background'].get('density_enforced',False),'source':p,**{k:d.get(k) for k in fields},'phase_inference_eligible':False})
t=read('results/physics48_torus.json');e=read('results/physics48_torus_excitations.json')
assert sum([5636,4*5454,4*5292])==48620
assert max(t['ground_vector_residuals'])<1e-10
assert max(o['full_commutator_max'] for o in t['operators'].values())<1e-10
assert abs(e['records'][0]['energies'][0]+7.342883075944)<1e-9
assert abs(e['records'][1]['energies'][0]+7.296046028855)<1e-9
out={'finite_states':rows,'matched_bond_changes':changes,'infinite_states':infinite,'accepted_gamma_fit_points':0,'gamma_fits':None,'reason_no_fit':'Unconverged entropy/energy, impure CGS winding sectors, unresolved bulk filling and MES; only two narrow circumferences.'}
(ROOT/'results/physics48_summary.json').write_text(json.dumps(out,indent=2)+'\n')
for r in rows:
 if r['cap']==512: print(r['family'],r['width'],r['owner_A_row_densities'],r['maximum_truncation_error'])
