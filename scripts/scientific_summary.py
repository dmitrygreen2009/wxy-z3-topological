"""Summarize measured states without promoting unconverged fits to topology."""
import glob,json,pathlib,itertools
import numpy as np
ln3=float(np.log(3))
finite=[]
for path in sorted(glob.glob('results/*_L*_w*_chi*.json')):
    r=json.load(open(path))
    if 'records' not in r or 'family' not in r or 'length' not in r:continue
    stages=r['records'];q=stages[-1]
    energies=q.get('sweep_energies',[])
    prior=[v for v in stages if v['cap']<q['cap']]
    finite.append(dict(file=path,family=r['family'],L=r['length'],w=r['width'],Ly=r['physical_circumference'],
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
    group=sorted(group,key=lambda r:r['Ly']);x=np.array([r['Ly'] for r in group]);y=np.array([r['S'] for r in group])
    alpha,b=np.polyfit(x,y,1);res=y-alpha*x-b
    alpha_fixed=float(x@(y+ln3)/(x@x));res_fixed=y-alpha_fixed*x+ln3
    pair_gammas=[]
    for i,j in itertools.combinations(range(len(x)),2):
        slope=(y[j]-y[i])/(x[j]-x[i]);pair_gammas.append(float(slope*x[i]-y[i]))
    fits.append(dict(family=family,ordering=ordering,L=L,widths=[r['w'] for r in group],chis=[r['chi'] for r in group],
        alpha=float(alpha),gamma_apparent=float(-b),gamma_minus_ln3=float(-b-ln3),
        residual_rms=float(np.sqrt(np.mean(res**2))),ln3_fixed_slope=alpha_fixed,
        ln3_fixed_residual_rms=float(np.sqrt(np.mean(res_fixed**2))),pairwise_gamma_range=[min(pair_gammas),max(pair_gammas)],
        interpretation='Descriptive finite-size fit. Bond, length, circumference and state-sector convergence must all hold before testing a topological constant.'))
infinite=[]
for path in sorted(glob.glob('results/infinite_*_chi*.json')):
    r=json.load(open(path))
    if 'spatial_entropy' not in r:continue
    r['file']=path;infinite.append(r)
output=dict(ln3=ln3,finite_states=finite,best_finite_states=best,entropy_fits=fits,infinite_states=infinite,
    status='Study in progress; transfer lengths describe the measured variational infinite states, and apparent gamma values are not converged thermodynamic estimates.')
pathlib.Path('results/scientific_summary.json').write_text(json.dumps(output,indent=2))
print('Finite states:',len(finite),'Infinite states:',len(infinite),'Descriptive fits:',len(fits))
for r in best:print(r['family'],r['L'],r['w'],r['ordering'],r['chi'],r['E'],r['S'])
