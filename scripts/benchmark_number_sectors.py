"""Ground-energy scan over every exact microscopic number sector.
This scan never estimates complete degeneracies from a small scalar Krylov list.
"""
import json,pathlib
from time import time
import numpy as np
from scipy.sparse.linalg import eigsh
from ed import hamiltonian
from provenance import provenance
models={'single_star':[[0,1,2]],'shared_edge_pair':[[0,1,2],[0,3,4]],'four_star_torus':[[0,1,2],[0,4,2],[3,4,5],[3,1,5]]}
results={}
for name,legs in models.items():
    n=3*len(legs)+max(e for es in legs for e in es)+1
    points=[]
    for nup in range(n+1):
        started=time();h,basis,herm=hamiltonian(legs,nup)
        audit=provenance(7113+nup,'Independent SciPy ED number-sector scan',dict(tol=1e-11,ncv=40,k=1,residual_tolerance=1e-8),'Complex random starting vector',[f'N_up={nup}'])
        if len(basis)<=4:
            values,vectors=np.linalg.eigh(h.toarray());e=float(values[0]);v=vectors[:,0]
        else:
            rng=np.random.default_rng(7113+nup)
            values,vectors=eigsh(h,k=1,which='SA',tol=1e-11,ncv=min(40,len(basis)),maxiter=10000,v0=rng.normal(size=len(basis))+1j*rng.normal(size=len(basis)))
            e=float(values[0]);v=vectors[:,0]
        residual=float(np.linalg.norm(h@v-e*v))
        assert herm<1e-13 and residual<1e-8
        points.append(dict(nup=nup,spins=n,dimension=len(basis),energy=e,residual=residual,hermiticity_error=herm,runtime_seconds=time()-started,audit=audit))
        results[name]=dict(points=points,minimum_energy=min(p['energy'] for p in points),complete_number_scan=len(points)==n+1,
            interpretation='Ground energies by number sector only; no multiplicities inferred from k=1.')
        pathlib.Path('results/benchmark_number_sectors.json').write_text(json.dumps(results,indent=2)+'\n')
        print(name,nup,len(basis),e,residual,flush=True)
    minimum=min(p['energy'] for p in points)
    results[name]['minimizing_number_sectors']=[p['nup'] for p in points if abs(p['energy']-minimum)<1e-8]
    assert max(abs(points[q]['energy']-points[n-q]['energy']) for q in range(n+1))<1e-8
pathlib.Path('results/benchmark_number_sectors.json').write_text(json.dumps(results,indent=2)+'\n')
