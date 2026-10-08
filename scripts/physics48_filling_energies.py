"""Inventory saved number-sector energies; never run a solver or infer minima from incomplete scans."""
import csv,json,pathlib,re,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=pathlib.Path(__file__).resolve().parents[1];records=[];sources={};unreadable=[]
def load(p):return json.loads(p.read_text())
for file in ['benchmark_number_sectors.json','small_cylinder_number_sectors.json']:
 p=ROOT/'results'/file;d=load(p);sources[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
 for geo,v in d.items():
  for x in v['points']:
   records.append(dict(method='ED',geometry=geo,nup=x['nup'],spins=x['spins'],filling=x['nup']/x['spins'],energy=x['energy'],energy_per_site=x['energy']/x['spins'],residual=x['residual'],cap=None,source=str(p.relative_to(ROOT)),record_path=f"{geo}.points[N={x['nup']}]",number_evidence='explicit fixed sector',branch='ED',maximum_truncation_error=None,last_sweep_energy_change=None))
# Only documented finite/infinite result payload schemas. Skip checkpoints,
# logs and diagnostic copies to avoid treating partial solver output as final.
measured={load(p)['source_result']:load(p)['nup'] for p in (ROOT/'results').glob('physics48_*_L*_w*_chi*.json') if not p.stem.endswith('_error')}
files=list((ROOT/'results').glob('*_star*.json'))+list((ROOT/'results').glob('*charge_sectors*.json'))+list((ROOT/'results'/'snapshots').glob('*number_audit*.json'))
files=sorted(set(files));seen=set()
def append_finite(x,p,context,path):
 merged={**context,**x};ns=merged.get('physical_spins',merged.get('spins'));n=merged.get('nup');e=x.get('energy');family=merged.get('family');L=merged.get('length');w=merged.get('width')
 if not all(isinstance(v,(int,float)) for v in [ns,e]) or not family or L is None or w is None:return
 source=str(p.relative_to(ROOT));evidence='explicit fixed sector'
 if n is None and source in measured:n=measured[source];evidence='saved checkpoint QN measurement'
 if n is None:
  audit=merged.get('audit',{});labels=audit.get('conserved_quantum_numbers',[])
  matches=[re.search(r'N_up\s*=\s*(\d+)',s) for s in labels];matches=[m for m in matches if m]
  if matches:n=int(matches[0].group(1));evidence='saved audit sector label'
 if n is None:
  # Historical neutral driver used ceil(N/2); label this weaker evidence.
  if 'delta_nup' not in merged and ('records' in context or 'records' in x):n=(int(ns)+1)//2;evidence='historical neutral-driver initialization; not remeasured'
  else:return
 sweeps=x.get('sweep_energies',[]);trunc=x.get('sweep_max_truncation_errors',x.get('truncation_errors',[]));cap=x.get('cap',merged.get('chi',merged.get('bond_dimension')))
 key=(source,path,n,cap,e)
 if key in seen:return
 seen.add(key);sources[source]=hashlib.sha256(p.read_bytes()).hexdigest()
 branch=merged.get('ordering','unknown')+':'+str(merged.get('seed',merged.get('audit',{}).get('random_seed','unknown')))
 records.append(dict(method='variational',geometry=f'{family}_L{L}_w{w}',nup=int(n),spins=int(ns),filling=n/ns,energy=e,energy_per_site=e/ns,residual=None,cap=cap,source=source,record_path=path,number_evidence=evidence,branch=branch,maximum_truncation_error=max(trunc) if trunc else None,last_sweep_energy_change=abs(sweeps[-1]-sweeps[-2]) if len(sweeps)>1 else None))
for p in files:
 try:d=load(p)
 except Exception as exc:unreadable.append({'source':str(p.relative_to(ROOT)),'error':str(exc)});continue
 if not isinstance(d,dict) or d.get('kind')=='infinite' or p.name.startswith('physics48_'):continue
 context=d.get('geometry',d)
 if not isinstance(context,dict):continue
 if 'records'in d:
  for i,r in enumerate(d['records']):
   if isinstance(r,dict):append_finite(r,p,context,f'records[{i}]')
 else:append_finite(d,p,context,'root')
minima=[]
for method in ['ED','variational']:
 for geo in sorted({r['geometry'] for r in records if r['method']==method}):
  rr=[r for r in records if r['method']==method and r['geometry']==geo];best=min(r['energy'] for r in rr);bestrows=[r for r in rr if abs(r['energy']-best)<1e-8]
  minima.append(dict(method=method,geometry=geo,minimum_recorded_energy=best,lowest_recorded_numbers=sorted(set(r['nup'] for r in bestrows)),lowest_recorded_fillings=sorted(set(r['filling'] for r in bestrows)),number_sectors_recorded=sorted(set(r['nup'] for r in rr)),complete_number_scan=method=='ED' and len(set(r['nup'] for r in rr))==rr[0]['spins']+1,qualification='Complete finite ED scan' if method=='ED' else 'Lowest stored variational upper bound only; branches and caps differ; no ground-sector certification'))
finite=[r for r in records if r['method']=='ED'];ph=[]
for geo in sorted(set(r['geometry'] for r in finite)):
 rr=[r for r in finite if r['geometry']==geo];by={r['nup']:r for r in rr};ns=rr[0]['spins'];err=max(abs(r['energy']-by[ns-r['nup']]['energy']) for r in rr);assert err<1e-8
 ph.append(dict(geometry=geo,particle_hole_energy_error=err,maximum_eigenvector_residual=max(r['residual'] for r in rr)))
resolved=[]
for name in ['winding_physical_sector_ed.json','sector_penalty_ed_audit.json']:
 p=ROOT/'results'/name;d=load(p);sources[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
 for r in d['records']:
  for sec in r['sectors']:
   resolved.append(dict(source=str(p.relative_to(ROOT)),geometry=f"{r['family']}_L{r['L']}_w{r['width']}",nup=r['nup'],spins=r['physical_spins'],charge=sec['charge'],energy=sec['energy'],energy_per_site=sec['energy']/r['physical_spins'],role='low-number regression fixture' if 'penalty' in name else 'physical winding-sector ground candidate'))
output=dict(symmetry_resolved_ED=resolved,records=records,minima=minima,exact_particle_hole_checks=ph,source_sha256=sources,unreadable_sources=unreadable,scope='All saved complete ED scans and finite result/charge-sector payloads matching documented energy schemas, including all recorded caps. No fabricated reflected variational measurements. Checkpoint iteration metadata and ad hoc energies in logs excluded. Infinite non-number-conserving expectation values remain in physics48_summary.json, not exact sector curves.',remaining_limitations='Missing saved results are not reconstructed from log text; finite variational minima cannot order unconverged sectors rigorously.')
output['analysis_git_commit']=__import__('subprocess').check_output(['git','rev-parse','HEAD'],text=True).strip()
(ROOT/'results/physics48_filling_energies.json').write_text(json.dumps(output,indent=2)+'\n')
with open(ROOT/'results/physics48_filling_energies.csv','w') as f:
 writer=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');writer.writeheader();writer.writerows(records)
fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for geo in sorted(set(r['geometry'] for r in finite)):
 rr=sorted([r for r in finite if r['geometry']==geo],key=lambda r:r['filling']);axes[0].plot([r['filling'] for r in rr],[r['energy_per_site'] for r in rr],'.-',label=geo)
for geo in sorted(set(r['geometry'] for r in records if r['method']=='variational')):
 rr=[r for r in records if r['method']=='variational' and r['geometry']==geo];best={}
 for r in rr:
  if r['nup'] not in best or r['energy']<best[r['nup']]['energy']:best[r['nup']]=r
 rr=sorted(best.values(),key=lambda r:r['filling'])
 if len(rr)<2:continue
 axes[1].scatter([r['filling'] for r in rr],[r['energy_per_site'] for r in rr],s=16,label=geo)
for ax,title in zip(axes,['Exact finite-sector minima','Saved variational upper bounds']):
 ax.set(xlabel='N_up / physical sites',ylabel='Energy / physical site',title=title);ax.axvline(.5,color='grey',lw=.6);ax.legend(fontsize=6);ax.grid(alpha=.2)
fig.savefig(ROOT/'results/physics48_filling_energies.png',dpi=160);plt.close(fig)
print('Energy records',len(records),'unreadable',len(unreadable));print(json.dumps(minima,indent=2))
