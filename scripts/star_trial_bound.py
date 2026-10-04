"""Independent ED validation of the isolated-A-star variational upper bound."""
import json,pathlib,time
import numpy as np
from ed import hamiltonian
from provenance import provenance
started=time.time()
hs,bs,herm=hamiltonian([[0,1,2]],3)
values,vectors=np.linalg.eigh(hs.toarray());ground=vectors[:,0]
amplitudes={int(b):v for b,v in zip(bs,ground)}
legs=[[0,1,2],[0,4,2],[3,4,5],[3,1,5]]
h,basis,herm=hamiltonian(legs,9)
psi=np.zeros(len(basis),complex)
A_indices=[[0,1,2,12,13,14],[6,7,8,15,16,17]]
for j,raw in enumerate(basis):
    b=int(raw)
    if (b>>3)&7 != 7 or (b>>9)&7 != 0:continue
    local=[sum(((b>>physical)&1)<<a for a,physical in enumerate(indices)) for indices in A_indices]
    psi[j]=amplitudes.get(local[0],0)*amplitudes.get(local[1],0)
assert abs(np.linalg.norm(psi)-1)<1e-12
energy=complex(np.vdot(psi,h@psi));target=2*values[0]
assert abs(energy-target)<1e-12
assert values[0]<-2.4030921
output=dict(physical_spins=18,nup=9,endpoint_leg_table_zero_based=legs,
    isolated_A_star_physical_indices_zero_based=A_indices,B_matter_configuration='vertex 2 all up; vertex 4 all down',
    isolated_star_energy=float(values[0]),trial_energy=energy.real,expectation_error=float(abs(energy-target)),
    norm=float(np.linalg.norm(psi)),isolated_star_residual=float(np.linalg.norm(hs@ground-values[0]*ground)),
    infinite_two_slice_bound_per_width=-4.8061842,
    derivation='Each of 2*width A stars is a disjoint six-spin fixed-Nup=3 state. All B matter spins are Z products, so all B exchange expectations vanish. Choose 3*width B matter up spins for half filling. Round the validated trial energy upward.',
    interpretation='A ground-state upper bound valid for either cylinder family; satisfying it is not convergence.',
    runtime_seconds=time.time()-started,
    audit=provenance(None,'Independent ED product-of-stars expectation',dict(expectation_tolerance=1e-12),'Exact 20-dimensional star eigenvector and deterministic B matter',['N_up=9']))
pathlib.Path('results/star_trial_bound_audit.json').write_text(json.dumps(output,indent=2)+'\n')
print('Star-product trial energy:',energy.real,'error:',abs(energy-target),'bound per infinite cell width:',-4.8061842)
