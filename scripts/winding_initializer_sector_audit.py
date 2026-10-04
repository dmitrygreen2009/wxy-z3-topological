"""Resolve additional exact CGS blocks behind the direct-QN initializer trap."""
import itertools,json,pathlib,time
import numpy as np
from provenance import provenance
start=time.time();family='zigzag';g=json.load(open('geometry/zigzag_L2_w1_star.json'));n=g['physical_spins'];nv=len(g['vertices'])
legs=[[None]*3 for _ in range(nv)]
for e in g['shared_gauge_edges']:
 for v in e['endpoints']:legs[v['vertex_id']-1][v['local_leg_i']-1]=e['edge_id']-1
basis=[sum(1<<j for j in c) for c in itertools.combinations(range(n),2)];lookup={b:j for j,b in enumerate(basis)}
local=json.load(open('results/matter_charge_basis_audit.json'));h=np.zeros((len(basis),len(basis)),complex)
for v,es in enumerate(legs):
 for i,e in enumerate(es):
  gauge=3*nv+e
  for term in local['rotated_exchange_operators'][i]['nonzero_matrix_units']:
   coef=complex(term['real'],term['imag'])
   for col,old in enumerate(basis):
    if (old>>gauge)&1 or ((old>>(3*v))&7)!=term['input']:continue
    new=(old&~(7<<(3*v)))|(term['out']<<(3*v))|(1<<gauge);row=lookup[new]
    h[row,col]-=coef;h[col,row]-=coef.conjugate()
cycles=json.load(open('geometry/cgs_cycles/zigzag_L2_w1_star.json'))['cycles'];charges=[]
for cycle in cycles:
 p=cycle['edge_exponents'];weights=[0]*n
 for v,es in enumerate(legs):
  c=p[es[0]]%3;q=(p[es[1]]-c)%3
  for a in range(3):weights[3*v+a]=(c+q*a)%3
 for e,pe in enumerate(p):weights[3*nv+e]=pe%3
 charges.append(np.array([sum(w*((b>>j)&1) for j,w in enumerate(weights))%3 for b in basis]))
records=[]
for q1 in range(3):
 chosen=np.flatnonzero((charges[0]==0)&(charges[1]==q1));vals=np.linalg.eigvalsh(h[np.ix_(chosen,chosen)])
 records.append(dict(requested_winding_q=0,additional_winding_q=q1,dimension=len(chosen),first_energies=vals[:8].tolist(),ground_energy=float(vals[0])))
print([(r['additional_winding_q'],r['ground_energy']) for r in records],flush=True)
out=dict(family=family,L=2,width=1,Nup=2,total_sz=2-n/2,physical_spins=n,records=records,
 audit=provenance(None,'Independent exact additional-CGS-block ED',dict(tolerance=1e-11),'Complete low-number blocks grouped by both microscopic winding operators',['N_up','both winding charges']),
 runtime_seconds=time.time()-start,
 interpretation='New initializer-trapping audit, not a repeated benchmark. Matching a restricted-block energy is not proof of the saved state charge; direct saved-state loop measurement remains required.')
pathlib.Path('results/winding_initializer_sector_audit.json').write_text(json.dumps(out,indent=2)+'\n')
