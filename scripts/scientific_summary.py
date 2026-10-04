"""Summarize measured states without promoting unconverged fits to topology."""
import glob,json,pathlib,itertools
import numpy as np
from finite_results import finite_groups
from fit_entropy import fit as entropy_fit
ln3=float(np.log(3))
finite=[]
for path,r,indexed_stages in finite_groups():
    stages=[q for k,q in indexed_stages];index,q=indexed_stages[-1]
    energies=q.get('sweep_energies',[])
    prior=[v for v in stages if v['cap']<q['cap']]
    finite.append(dict(file=path,record_index=index,family=r['family'],L=r['length'],w=r['width'],Ly=r['physical_circumference'],
        ordering=r.get('ordering','axial'),nup=r.get('nup',(r['spins']+1)//2),ground_state_scope='Fixed number sector; global cylinder minimum not certified',spins=r['spins'],chi=q['cap'],actual_chi=q['maxlinkdim'],E=q['energy'],S=q['entropy'],
        last_sweep_energy_change=energies[-1]-energies[-2] if len(energies)>1 else None,
        last_truncation_error=q.get('sweep_max_truncation_errors',[None])[-1],
        entropy_refinement_change=q.get('entropy_change_refinement'),
        entropy_change_from_lower_chi=q['entropy']-prior[-1]['entropy'] if prior else None,
        energy_change_from_lower_chi=q['energy']-prior[-1]['energy'] if prior else None))
# One record per geometry/order: select lowest energy across sectors and seeds.
best=[]
keys=sorted(set((r['family'],r['L'],r['w'],r['ordering']) for r in finite))
for key in keys:
    candidates=[r for r in finite if (r['family'],r['L'],r['w'],r['ordering'])==key]
    best.append(sorted(candidates,key=lambda r:(r['E'],-r['chi']))[0])
fits=[]
for family,ordering,L in sorted(set((r['family'],r['ordering'],r['L']) for r in best)):
    group=[r for r in best if r['family']==family and r['ordering']==ordering and r['L']==L]
    if len(group)<3:continue
    group=sorted(group,key=lambda r:r['Ly'])
    points=[dict(circumference=r['Ly'],entropy=r['S'],width=r['w']) for r in group]
    estimate=entropy_fit(points)
    estimate.update(family=family,ordering=ordering,L=L,widths=[r['w'] for r in group],chis=[r['chi'] for r in group],
        raw_points=group,smallest_circumference_removed=entropy_fit(points[1:]),
        eligible_for_topological_inference=False,
        interpretation='Descriptive lowest-energy envelope across recorded caps and number/flux branches. No common MES or global-ground selection is certified; no thermodynamic gamma confidence interval.')
    fits.append(estimate)

bond_audit={r["source_result"]:r for r in json.load(open("results/infinite_bond_dimension_audit.json"))["records"]}
infinite=[]
for path in sorted(glob.glob('results/infinite_*_chi*.json')):
    r=json.load(open(path))
    if 'spatial_entropy' not in r or r.get('validation_fixture',False):continue
    r={k:v for k,v in r.items() if not any(term in k.lower() for term in ('triplet','charge_three'))}
    r['excluded_observable_policy']='Matter-triplet measurements are archival only'
    r['file']=path
    r['historical_reported_chi']=r['chi']
    r['audited_maximum_bond_dimension']=bond_audit.get(path,{}).get('maximum_bond_dimension',r.get('bond_dimension_details',{}).get('maximum_bond_dimension'))
    r['bond_dimension_audit_status']=bond_audit.get(path,{}).get('status','Not yet independently audited')
    infinite.append(r)
output=dict(ln3=ln3,finite_states=finite,best_finite_states=best,entropy_fits=fits,infinite_states=infinite,
    scientific_strategy_reaudit='results/scientific_strategy_reaudit.md',conditional_filling_constraints=json.load(open('results/primitive_filling_audit.json')),phase_claims=dict(Z3_topological_order='not established',Z3_topological_order_excluded='not established',two_dimensional_gaplessness='not established',physical_U1_order='not established'),status='Study in progress; transfer lengths describe the measured variational infinite states, and apparent gamma values are not converged thermodynamic estimates.')
pathlib.Path('results/scientific_summary.json').write_text(json.dumps(output,indent=2))
print('Finite states:',len(finite),'Infinite states:',len(infinite),'Descriptive fits:',len(fits))
for r in best:print(r['family'],r['L'],r['w'],r['ordering'],r['chi'],r['E'],r['S'])
