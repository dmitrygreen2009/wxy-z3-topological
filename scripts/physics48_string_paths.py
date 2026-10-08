"""New same-endpoint path test only; reuse cached torus states, no H rebuild or solver."""
import json,pathlib,time
import numpy as np
start=time.time();cache=np.load('results/physics48_torus_joint_vectors.npz');basis=cache['basis'];G=cache['vectors'];labels=json.load(open('results/physics48_torus.json'))['joint_cycle_states'];lookup={int(b):i for i,b in enumerate(basis)}
legs=[[0,1,2],[0,4,2],[3,4,5],[3,1,5]];omega=np.exp(2j*np.pi/3)
def apply(p,V):
 mapped=basis.copy();phase=np.zeros(len(basis),dtype=int)
 for v,es in enumerate(legs):
  c=p[es[0]];q=(p[es[1]]-c)%3;mapped &= ~(7<<(3*v))
  for a in range(3):
   bit=(basis>>(3*v+a))&1;mapped|=bit<<(3*v+(a+q)%3);phase+=c*bit
 for e,x in enumerate(p):phase+=x*((basis>>(12+e))&1)
 output=np.empty_like(V);output[[lookup[int(b)] for b in mapped]]=omega**(phase[:,None]%3)*V
 return output
short=[0,1,0,0,0,0];long=[1,0,0,1,2,0];closed=[1,2,0,1,2,0]
identity=np.linalg.norm(apply(long,G)-apply(short,apply(closed,G)))
charge_errors=[float(np.linalg.norm(apply(long,G[:,[k]])-omega**row['X_q']*apply(short,G[:,[k]]))) for k,row in enumerate(labels)]
assert identity<1e-12 and max(charge_errors)<1e-9
r=dict(short_path_edges=[1],long_path_edges=[0,4,3],same_endpoint_vertices=[0,3],difference_closed_cycle='X',full_monomial_relation='U_long=U_short U_X',ground_vector_relation_error=float(identity),ground_charge_phase_errors=charge_errors,interpretation='Path-independent up to an exact cycle eigenvalue for fixed endpoints and a definite cycle-sector state. This is an exact symmetry identity, not an R-dependent defect energy or a topological excitation identification.',runtime_seconds=time.time()-start)
pathlib.Path('results/physics48_string_paths.json').write_text(json.dumps(r,indent=2)+'\n');print(r)
