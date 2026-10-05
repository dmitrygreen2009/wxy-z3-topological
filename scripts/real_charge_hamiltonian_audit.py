"""Independent full eight-state proof of the real microscopic exchange form."""
import json, pathlib, time
import numpy as np
from ed import W
from provenance import provenance

started=time.time()
B=np.zeros((8,8),complex);B[0,0]=B[7,7]=1
for a in range(3):
    for b in range(3):
        B[1<<b,1<<a]=np.exp(-2j*np.pi*a*b/3)/np.sqrt(3)
        B[7^(1<<b),7^(1<<a)]=np.exp(2j*np.pi*a*b/3)/np.sqrt(3)
lower=[];z=[]
for a in range(3):
    op=np.zeros((8,8))
    for old in range(8):
        if (old>>a)&1:op[old^(1<<a),old]=1
    lower.append(op)
    z.append(np.diag([((old>>a)&1)-.5 for old in range(8)]))
alpha=(1-1/np.sqrt(3))/2;beta=2*(1+1/np.sqrt(3));delta=2/np.sqrt(3)
records=[]
for i in range(3):
    j,k=[a for a in range(3) if a!=i]
    physical=sum(W[a,i]*lower[a] for a in range(3))
    exact=B.conj().T@physical@B
    closed=alpha*lower[i]+beta*lower[i]@z[j]@z[k]+delta*lower[i].T@lower[j]@lower[k]
    error=float(np.linalg.norm(exact-closed));assert error<1e-13
    records.append(dict(local_leg_i=i+1,full_eight_state_operator_error=error,
        maximum_imaginary_roundoff=float(np.max(np.abs(exact.imag)))))
out=dict(formula='A_i = alpha S_i^- + beta S_i^- Sz_j Sz_k + delta S_i^+ S_j^- S_k^-; j,k are the other matter modes',
    coefficients=dict(alpha=alpha,beta=beta,delta=delta),records=records,
    matter_states_retained=8,product_terms_per_star_including_hc=18,
    runtime_seconds=time.time()-started,
    audit=provenance(None,'Independent full eight-state matrix identity',dict(tolerance=1e-13),'All local matter kets',['matter number flow','CGS charge flow']),
    interpretation='Exact real spin Hamiltonian in a local number-preserving basis. The three-mode term is required by the Hamiltonian basis change; it is not an added observable, interaction, clock approximation, or phase diagnostic. Production MPO equivalence tests are a separate gate.')
pathlib.Path('results/real_charge_hamiltonian_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print('Exact real exchange identity validated on all eight matter states for all three legs')
