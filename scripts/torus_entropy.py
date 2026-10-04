"""Full-state spatial entropy within the resolved microscopic torus groundspace."""
import json,pathlib
import numpy as np
from scipy.optimize import minimize

d=np.load('results/torus_groundspace.npz');basis=d['basis'];vectors=d['ground_vectors']
# Matter of A(0), B(0), their leg-1/3 gauge edges, and the leg-2 edge
# originating at A(1). The last shared edge is assigned to this side of cut.
left=[0,1,2,3,4,5,12,14,16];right=[i for i in range(18) if i not in left]
li=sum(((basis>>s)&1)<<k for k,s in enumerate(left))
ri=sum(((basis>>s)&1)<<k for k,s in enumerate(right))
matrices=np.zeros((4,512,512),dtype=complex)
for k in range(4):matrices[k,li,ri]=vectors[:,k]
sectors=[np.array([i for i in range(512) if i.bit_count()==q]) for q in range(10)]
blocks=[matrices[:,l,:][:,:,r] for l,r in zip(sectors,reversed(sectors))]

def probabilities(c):
    c=np.asarray(c,dtype=complex);c/=np.linalg.norm(c)
    p=np.concatenate([np.linalg.svd(np.einsum('k,kij->ij',c,b),compute_uv=False)**2 for b in blocks])
    return np.sort(p/p.sum())[::-1]
def entropy(c):
    p=probabilities(c);p=p[p>1e-15];return float(-sum(p*np.log(p)))
def coeffs(x):return x[:4]+1j*x[4:]

rng=np.random.default_rng(7103)
samples=[entropy(np.eye(4)[k]) for k in range(4)]
random_samples=[entropy(rng.normal(size=4)+1j*rng.normal(size=4)) for _ in range(12)]
x=rng.normal(size=8)
opt=minimize(lambda x:entropy(coeffs(x)),x,method='Powell',options={'maxiter':40,'maxfev':3000,'xtol':1e-6,'ftol':1e-9})
out={'left_physical_sites_zero_based':left,'right_physical_sites_zero_based':right,
     'ed_basis_vector_entropies':samples,'random_groundspace_entropies':random_samples,
     'optimized_entropy':float(opt.fun),'optimization_success':bool(opt.success),'function_evaluations':opt.nfev,
     'optimized_schmidt_probabilities':probabilities(coeffs(opt.x)).tolist(),
     'interpretation':'Tiny-torus full-state entropy and groundspace dependence only; no topological inference or certified global MES minimum.'}
pathlib.Path('results/torus_entropy.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
