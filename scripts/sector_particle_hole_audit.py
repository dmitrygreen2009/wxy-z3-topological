"""New antiunitary covariance check within the winding-sector construction.
No ground-state diagonalization is repeated.
"""
import json,pathlib,time
import numpy as np
from ed import hamiltonian
from provenance import provenance
started=time.time();records=[];omega=np.exp(2j*np.pi/3)
def loop_matrix(legs,p,basis):
    nv=len(legs);lookup={int(b):i for i,b in enumerate(basis)};U=np.zeros((len(basis),len(basis)),complex)
    for col,old in enumerate(basis):
        old=int(old);new=old;exponent=0
        for v,es in enumerate(legs):
            c=int(p[es[0]])%3;q=int(p[es[1]]-c)%3;new &= ~(7<<(3*v))
            for a in range(3):
                bit=(old>>(3*v+a))&1;new|=bit<<(3*v+(a+q)%3);exponent+=c*bit
        exponent+=sum(int(p[e])*((old>>(3*nv+e))&1) for e in range(len(p)))
        U[lookup[new],col]=omega**(exponent%3)
    return U
for saved in json.load(open('results/sector_penalty_ed_audit.json'))['records']:
    g=json.load(open(saved['geometry']));nv=len(g['vertices']);n=g['physical_spins'];legs=[[None]*3 for _ in range(nv)]
    for edge in g['shared_gauge_edges']:
        for end in edge['endpoints']:legs[end['vertex_id']-1][end['local_leg_i']-1]=edge['edge_id']-1
    N=saved['nup'];p=saved['edge_exponents'];h,b,_=hamiltonian(legs,N);hh,bh,_=hamiltonian(legs,n-N)
    lookup={int(v):j for j,v in enumerate(bh)};X=np.zeros((len(bh),len(b)),complex)
    for j,v in enumerate(b):X[lookup[((1<<n)-1)^int(v)],j]=1
    U=loop_matrix(legs,p,b);Uh=loop_matrix(legs,p,bh);phase=omega**(sum(p)%3)
    Herror=np.linalg.norm(hh@X-X@h.conjugate());Uerror=np.linalg.norm(Uh@X-phase*X@U.conjugate())
    assert max(Herror,Uerror)<1e-11
    records.append(dict(family=saved['family'],geometry=saved['geometry'],spins=n,number_sectors=[N,n-N],
        phase_exponent=sum(p)%3,particle_hole_hamiltonian_error=float(Herror),particle_hole_loop_error=float(Uerror),
        charge_mapping={str(q):(sum(p)-q)%3 for q in range(3)},diagonalization_repeated=False))
r=dict(records=records,runtime_seconds=time.time()-started,
    audit=provenance(None,'Independent microscopic antiunitary matrix covariance',dict(tolerance=1e-11),'Complete paired low-number blocks',['N_up']),
    conclusion='Particle-hole maps (N,q) to (Ns-N,sum(p)-q mod 3); it equates paired energies, but does not pin a minimum to half filling.')
pathlib.Path('results/sector_particle_hole_audit.json').write_text(json.dumps(r,indent=2)+'\n')
print('Both microscopic winding-sector particle-hole identities passed; no diagonalizations repeated')
