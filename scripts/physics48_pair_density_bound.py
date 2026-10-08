"""One additional eleven-spin operator inequality, reusing the star multiplier.
No lambda optimization, no ground benchmark repetition, and no MPS optimization.
"""
import json,pathlib,time
import numpy as np
from ed import hamiltonian
start=time.time();b=json.load(open('results/physics48_density_bound.json'));trial=json.load(open('results/physics48_tiled_density_bound.json'));lam=b['cluster_multiplier'];records=[]
for n in range(12):
 H,basis,herm=hamiltonian([[0,1,2],[0,3,4]],n)
 Q=np.array([(int(x)&63).bit_count()+((int(x)>>6)&1)+.5*(int(x)>>7).bit_count() for x in basis])
 A=H.toarray()-lam*np.diag(Q);e,V=np.linalg.eigh(A);rec=np.linalg.norm(A-(V*e)@V.conj().T);orth=np.linalg.norm(V.conj().T@V-np.eye(len(basis)));margin=rec+np.linalg.norm(A)*orth+1e-10
 records.append(dict(auxiliary_cluster_nup=n,dimension=len(basis),lowest_shifted_eigenvalue=float(e[0]),conservative_epsilon=float(e[0]-margin),reconstruction_error=float(rec),orthogonality_error=float(orth),safety_margin=float(margin)))
assert sum(r['dimension'] for r in records)==2048
eps=min(r['conservative_epsilon'] for r in records);lo=(trial['upward_rounded_trial_upper']-eps)/(9*lam)
assert lo>=trial['same_local_cluster_density_window'][0]-1e-9
r=dict(lambda_reused=lam,epsilon_pair_conservative=eps,density_window=[lo,1-lo],records=records,proof='H=sum_cells(h_A+h_B)>=M*epsilon_pair(lambda)+lambda*N_up; each pair has one internal physical gauge, four external physical gauges, and six matter spins. Q_pair=Nm+N_internal+.5*N_external. Shared external gauges belong to two overlapping pair operators and count once in total Q. Full-plane primitive leg0 perfect matching only; not a parallel-edge tiny quotient.',qualification='Exact operator inequality with a residual-checked floating-point endpoint, explicit margins, and no interval-arithmetic claim. Shifted-operator N labels are not the preferred filling of original H.',original_H_unchanged=True,runtime_seconds=time.time()-start)
pathlib.Path('results/physics48_pair_density_bound.json').write_text(json.dumps(r,indent=2)+'\n');print('Pair-based density window',r['density_window'],'epsilon',eps,'seconds',r['runtime_seconds'])
