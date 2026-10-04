"""Verify all nine local microscopic CGS monomial transformations, without ED."""
import json,pathlib,time
import numpy as np
from ed import W
from provenance import provenance
started=time.time();H=np.zeros((64,64),complex)
for b in range(64):
 for a in range(3):
  for i in range(3):
   m=(b>>a)&1;g=(b>>(3+i))&1
   if m!=g:H[b^(1<<a)^(1<<(3+i)),b]=-W[a,i] if m else -W[a,i].conjugate()
records=[]
for c in range(3):
 for q in range(3):
  p=[(c+q*i)%3 for i in range(3)];U=np.zeros((64,64),complex)
  for old in range(64):
   matter=sum(((old>>a)&1)<<((a+q)%3) for a in range(3));new=matter+(old&56)
   exponent=c*(old&7).bit_count()+sum(p[i]*((old>>(3+i))&1) for i in range(3))
   U[new,old]=np.exp(2j*np.pi*exponent/3)
  comm=float(np.linalg.norm(H@U-U@H));cube=float(np.linalg.norm(U@U@U-np.eye(64)))
  unitary=float(np.linalg.norm(U.conj().T@U-np.eye(64)))
  assert max(comm,cube,unitary)<1e-13
  records.append(dict(c=c,q=q,gauge_exponents=p,commutator_error=comm,cube_error=cube,unitarity_error=unitary))
output=dict(records=records,all_nine_transformations_passed=True,physical_spins=6,bit_convention='1 is up',
 runtime_seconds=time.time()-started,interpretation='Exact local symmetry identities, not a ground-state calculation or a topology claim.',
 audit=provenance(None,'Independent microscopic CGS matrix identities',dict(tolerance=1e-13),'All 64 physical kets, all nine affine Z3 endpoint exponent patterns',[]))
pathlib.Path('results/cgs_general_audit.json').write_text(json.dumps(output,indent=2)+'\n')
print('Nine local CGS transformations passed; maximum commutator error:',max(r['commutator_error'] for r in records))
