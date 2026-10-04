"""Check global small-cluster minima and charge gaps independently of DMRG."""
import json,pathlib
import numpy as np
from scipy.sparse.linalg import eigsh
from ed import hamiltonian
out={}
for name,legs in [('star',[[0,1,2]]),('pair',[[0,1,2],[0,3,4]]),('torus',[[0,1,2],[0,4,2],[3,4,5],[3,1,5]])]:
    n=3*len(legs)+max(sum(legs,[]))+1
    energies=[]
    for nup in range(n//2+1):
        h,b,err=hamiltonian(legs,nup)
        e=float(np.linalg.eigvalsh(h.toarray())[0]) if len(b)<1000 else float(eigsh(h,k=1,which='SA',tol=1e-11,maxiter=10000,return_eigenvectors=False)[0])
        energies.append(e);print(name,nup,e,flush=True)
    assert np.argmin(energies)==n//2
    out[name]={'spins':n,'nup_sectors':list(range(n//2+1)),'ground_energies':energies,
        'charge_gap_at_half_filling':energies[-2]-energies[-1],
        'upper_sectors':'Related by spin flip plus complex conjugation'}
    pathlib.Path('results/sectors.json').write_text(json.dumps(out,indent=2))
