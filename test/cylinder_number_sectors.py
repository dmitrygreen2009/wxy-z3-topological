"""Regression: the microscopic tiny armchair cylinder minimizes off half filling."""
import json,pathlib,sys,time
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from scipy.sparse.linalg import eigsh
from ed import hamiltonian
from provenance import provenance
geometry=json.loads(pathlib.Path('geometry/armchair_L2_w1_star.json').read_text())
legs=[[int(e)-1 for e in row] for row in geometry['endpoint_leg_table']]
references={8:-8.037741645057586,10:-7.994302046793339,12:-8.037741645057586}
records=[]
for nup,reference in references.items():
    started=time.time();seed=7118+nup;rng=np.random.default_rng(seed)
    H,basis,hermiticity=hamiltonian(legs,nup)
    e,V=eigsh(H,k=1,which='SA',tol=1e-11,ncv=40,maxiter=10000,v0=rng.normal(size=len(basis))+1j*rng.normal(size=len(basis)))
    energy=float(e[0]);residual=float(np.linalg.norm(H@V[:,0]-energy*V[:,0]))
    assert hermiticity<1e-14 and residual<1e-8 and abs(energy-reference)<1e-8
    records.append(dict(nup=nup,energy=energy,residual=residual,hermiticity_error=hermiticity,dimension=len(basis),runtime_seconds=time.time()-started,
        audit=provenance(seed,'Independent microscopic ED sector regression',dict(k=1,tol=1e-11,ncv=40,maxiter=10000),'Complex Gaussian start',[f'N_up={nup}'])))
assert records[0]['energy']<records[1]['energy']-0.04
assert abs(records[0]['energy']-records[2]['energy'])<1e-8
pathlib.Path('results/cylinder_number_sector_regression.json').write_text(json.dumps(dict(records=records,
    interpretation='Sector-energy regression only; complete global scan and independent block multiplicities are separate audits.'),indent=2))
print('Passed armchair off-half-filling sector-energy and particle-hole regression.')
