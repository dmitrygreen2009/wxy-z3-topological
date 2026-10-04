"""Check symmetry restriction and scalar ARPACK undercounting of degeneracy."""
import json,pathlib
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from ed import hamiltonian
H,basis,_=hamiltonian([[0,1,2],[0,4,2],[3,4,5],[3,1,5]],9)
ref=np.array([-7.364534839426]*2+[-7.342883075944,-7.337288472473]+[-7.296046028855]*2+[-7.289970027036]*2)
lookup={int(b):i for i,b in enumerate(basis)}
# Exact unitary: flip all spins, then exchange local matter species 2 and 3.
mapped=((1<<18)-1)^basis
for v in range(4):
    a=3*v+1;b=3*v+2
    difference=((mapped>>a)^(mapped>>b))&1
    mapped ^= (difference<<a)|(difference<<b)
partner=np.array([lookup[int(b)] for b in mapped])
assert np.all(partner[partner]==np.arange(len(basis)))
U=coo_matrix((np.ones(len(basis)),(partner,np.arange(len(basis)))),shape=H.shape).tocsr()
comm=H@U-U@H;assert max(abs(comm.data),default=0)<1e-13
representatives=np.where(np.arange(len(basis))<partner)[0]
assert len(representatives)==24310
parities=[]
for sign in [1,-1]:
    columns=np.tile(np.arange(len(representatives)),2)
    rows=np.concatenate([representatives,partner[representatives]])
    P=coo_matrix((np.concatenate([np.ones(len(representatives)),sign*np.ones(len(representatives))])/np.sqrt(2),(rows,columns)),shape=(len(basis),len(representatives))).tocsr()
    hp=P.conj().T@H@P
    e=np.sort(eigsh(hp,k=12,which='SA',ncv=80,tol=1e-12,return_eigenvectors=False))
    parities.append({'parity':sign,'dimension':hp.shape[0],'energies':e.tolist(),'matches_supplied_first8':bool(max(abs(e[:8]-ref))<1e-9)})
    print('parity',parities[-1],flush=True)

scalar=[]
for tol,ncv,seed in [(t,c,s) for t in [1e-6,1e-8,1e-10,1e-12] for c in [20,40] for s in [0,1]]:
    rng=np.random.default_rng(seed);start=rng.normal(size=len(basis))+1j*rng.normal(size=len(basis))
    e,V=eigsh(H,k=8,which='SA',ncv=ncv,tol=tol,v0=start,maxiter=10000)
    ix=np.argsort(e);e=e[ix];V=V[:,ix]
    residual=np.linalg.norm(H@V-V*e,axis=0)
    r={'tol':tol,'ncv':ncv,'seed':seed,'energies':e.tolist(),'residual_max':float(max(residual)),
        'returned_ground_vectors':int(sum(abs(e-e[0])<1e-9)),
        'matches_supplied_first8':bool(max(abs(e-ref))<1e-9)}
    scalar.append(r);print('scalar',r,flush=True)
    pathlib.Path('results/torus_solver_audit.json').write_text(json.dumps({'parity_blocks':parities,'scalar_arpack':scalar},indent=2))
