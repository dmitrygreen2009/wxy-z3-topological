import json,pathlib
import numpy as np
from scipy.sparse.linalg import eigsh
from ed import hamiltonian
geometry=json.load(open('results/small_cylinder_incidence.json'));out={}
for family,g in geometry.items():
    n=g['spins'];nup=(n+1)//2
    H,basis,herm=hamiltonian(g['legs_zero_based'],nup)
    e,V=eigsh(H,k=4,which='SA',tol=1e-11,ncv=60,maxiter=10000)
    p=np.argsort(e);e=e[p];V=V[:,p]
    residual=np.linalg.norm(H@V-V*e,axis=0)
    assert max(residual)<1e-8
    left=g['left_physical_sites_zero_based'];right=[i for i in range(n) if i not in left]
    li=sum(((basis>>s)&1)<<k for k,s in enumerate(left));ri=sum(((basis>>s)&1)<<k for k,s in enumerate(right))
    M=np.zeros((2**len(left),2**len(right)),complex);M[li,ri]=V[:,0]
    schmidt=[]
    for q in range(len(left)+1):
        ls=[i for i in range(2**len(left)) if i.bit_count()==q]
        rs=[i for i in range(2**len(right)) if i.bit_count()==nup-q]
        if ls and rs:schmidt.extend(np.linalg.svd(M[np.ix_(ls,rs)],compute_uv=False)**2)
    schmidt=np.array(schmidt);schmidt/=sum(schmidt);nz=schmidt[schmidt>1e-15]
    out[family]={'spins':n,'nup':nup,'dimension':len(basis),'energies':e.tolist(),
        'residuals':residual.tolist(),'hermiticity_error':herm,'full_state_entropy_one_ED_eigenvector':float(-sum(nz*np.log(nz))),
        'entropy_interpretation':'State choice can differ from DMRG within a degenerate groundspace; energy is the independent regression target.'}
    pathlib.Path('results/small_cylinders_ed.json').write_text(json.dumps(out,indent=2));print(family,out[family],flush=True)
