"""Verify exact microscopic Z3 cycle symmetries, including local matter action.

For edge exponents p with sum incident p=0 mod 3, set c_v=p_leg1,
q_v=p_leg2-p_leg1. Rotate all local matter by c_v, cyclically permute
their labels by q_v, and rotate gauge spins by their p_e. Spin-up-number
rotations use omega**(p*n_up), fixing an irrelevant global phase so U^3=I.
"""
import json,pathlib
import numpy as np
from scipy.sparse import coo_matrix,eye
from ed import hamiltonian

legs=[[0,1,2],[0,4,2],[3,4,5],[3,1,5]]
h,basis,_=hamiltonian(legs,9);lookup={int(b):i for i,b in enumerate(basis)}
d=np.load('results/torus_groundspace.npz');G=d['ground_vectors']
assert np.array_equal(d['basis'],basis)
omega=np.exp(2j*np.pi/3)
cycles=[[1,0,2,0,0,0],[0,1,2,0,1,2],[0,0,0,1,0,2]]
out=[];projected=[]
for p in cycles:
    assert all(sum(p[e] for e in es)%3==0 for es in legs)
    mapped=basis.copy();exponent=np.zeros(len(basis),dtype=np.int64)
    for v,es in enumerate(legs):
        c=p[es[0]];q=(p[es[1]]-c)%3
        mapped &= ~(7<<(3*v))
        for a in range(3):
            bit=(basis>>(3*v+a))&1
            mapped |= bit<<(3*v+(a+q)%3)
            exponent += c*bit
    for e,pe in enumerate(p):exponent += pe*((basis>>(12+e))&1)
    rows=np.array([lookup[int(b)] for b in mapped]);phase=omega**(exponent%3)
    U=coo_matrix((phase,(rows,np.arange(len(basis)))),shape=h.shape).tocsr()
    def err(x):return float(max(abs(x.data),default=0))
    commute=err(h@U-U@h);order3=err(U@U@U-eye(len(basis)))
    assert commute<1e-13 and order3<1e-13
    g=G.conj().T@(U@G);projected.append(g)
    leak=float(np.linalg.norm(U@G-G@g));assert leak<1e-8
    out.append({'edge_exponents':p,'commutator_max_abs':commute,'U_cubed_error':order3,'groundspace_leakage':leak})

joint=projected[0]+0.317*projected[1]+0.137j*projected[2]
vals,V=np.linalg.eig(joint);V/=np.linalg.norm(V,axis=0)
charges=[]
for k in range(4):
    c=[]
    for u in projected:
        z=np.vdot(V[:,k],u@V[:,k]);nearest=int(np.argmin(abs(z-omega**np.arange(3))))
        assert abs(z-omega**nearest)<1e-8;c.append(nearest)
    charges.append(c)
result={'generators':out,'resolved_groundspace_cycle_charges':charges,
    'interpretation':'Exact microscopic commuting Z3 cycle symmetries on the tiny multigraph; this alone is not evidence of deconfinement or D(Z3) topological order. No effective clock/gauge Hamiltonian was used.'}
pathlib.Path('results/loop_symmetries.json').write_text(json.dumps(result,indent=2));print(result,flush=True)
