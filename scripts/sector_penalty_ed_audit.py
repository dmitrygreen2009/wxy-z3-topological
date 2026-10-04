"""Independent microscopic ED audit of winding projectors and exact penalties.
Low-number blocks validate algebra, not the ground filling or two-dimensional phase.
"""
import json,pathlib,time
import numpy as np
from ed import hamiltonian,W
from provenance import provenance
started=time.time();records=[]
for family in ('zigzag','armchair'):
    geometry=f'geometry/{family}_L2_w1_star.json'
    g=json.load(open(geometry));nv=len(g['vertices']);legs=[[None]*3 for _ in range(nv)]
    for edge in g['shared_gauge_edges']:
        for endpoint in edge['endpoints']:
            legs[endpoint['vertex_id']-1][endpoint['local_leg_i']-1]=edge['edge_id']-1
    assert all(None not in row for row in legs)
    cycle=json.load(open(f'geometry/cgs_cycles/{family}_L2_w1_star.json'))['cycles'][0]
    p=np.array(cycle['edge_exponents']);assert cycle['winding_number']==1
    h,basis,_=hamiltonian(legs,2);h=h.toarray();lookup={int(b):i for i,b in enumerate(basis)}
    U=np.zeros_like(h);omega=np.exp(2j*np.pi/3)
    for col,old in enumerate(basis):
        old=int(old);new=old;exponent=0
        for v,es in enumerate(legs):
            c=int(p[es[0]])%3;q=int(p[es[1]]-c)%3
            assert [(c+q*i)%3 for i in range(3)]==[int(p[e])%3 for e in es]
            mask=7<<(3*v);new &= ~mask
            for a in range(3):
                bit=(old>>(3*v+a))&1
                new |= bit<<(3*v+(a+q)%3)
                exponent+=c*bit
        exponent+=sum(int(p[e])*((old>>(3*nv+e))&1) for e in range(len(p)))
        U[lookup[new],col]=omega**(exponent%3)
    identity=np.eye(len(basis));comm=np.linalg.norm(h@U-U@h)
    cube=np.linalg.norm(U@U@U-identity);unitarity=np.linalg.norm(U.conj().T@U-identity)
    assert max(comm,cube,unitarity)<1e-11
    bound=nv*np.abs(W).sum();strength=2*bound+1;sectors=[]
    for charge in range(3):
        z=omega**charge;P=(identity+z.conjugate()*U+z.conjugate()**2*(U@U))/3
        assert np.linalg.norm(P-P.conj().T)<1e-11
        assert np.linalg.norm(P@P-P)<1e-11
        values,vectors=np.linalg.eigh(P);V=vectors[:,values>1-1e-10]
        assert V.shape[1]>0
        target,Vsmall=np.linalg.eigh(V.conj().T@h@V);state=V@Vsmall[:,0]
        penalized=h+strength*(identity-P);ep,vp=np.linalg.eigh(penalized)
        normres=np.linalg.norm(h@state-target[0]*state);mean=np.vdot(state,U@state)
        assert normres<1e-10
        assert abs(ep[0]-target[0])<1e-10
        assert abs(np.vdot(vp[:,0],U@vp[:,0])-z)<1e-10
        assert np.linalg.norm((penalized-h)@V)<1e-10
        sectors.append(dict(charge=charge,dimension=V.shape[1],energy=float(target[0]),penalized_energy=float(ep[0]),
            eigenvector_residual=float(normres),loop_mean_real=float(mean.real),loop_mean_imag=float(mean.imag),
            unitary_loop_variance=float(1-abs(mean)**2),penalty_within_sector_error=float(np.linalg.norm((penalized-h)@V))))
    records.append(dict(family=family,L=2,width=1,physical_spins=g['physical_spins'],nup=2,
        total_sz=2-g['physical_spins']/2,filling=2/g['physical_spins'],number_sector='Fixed fixture, not variational filling selection',
        geometry=geometry,edge_exponents=p.tolist(),commutator_error=float(comm),cube_error=float(cube),unitarity_error=float(unitarity),
        norm_bound=float(bound),penalty_strength=float(strength),sectors=sectors))
    print(family,'all three sectors validated',flush=True)
audit=provenance(7251,'Independent NumPy full ED of low-number microscopic blocks',dict(tolerance=1e-10),'Complete N_up=2 bases',['N_up'])
pathlib.Path('results/sector_penalty_ed_audit.json').write_text(json.dumps(dict(records=records,audit=audit,runtime_seconds=time.time()-started,
    interpretation='Exact loop/penalty validation only; ground filling, optimized MPS convergence and topological flux/MES identification are separate gates.'),indent=2)+'\n')
