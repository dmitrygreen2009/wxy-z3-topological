"""Independent full-number-block audit of the exact rotated spin Hamiltonian."""
import json,pathlib,time
import numpy as np
from ed import hamiltonian
from provenance import provenance
start=time.time();local=json.load(open('results/matter_charge_basis_audit.json'))
B=np.array(local['basis_real'])+1j*np.array(local['basis_imag']);records=[]
for saved in json.load(open('results/sector_penalty_ed_audit.json'))['records']:
    g=json.load(open(saved['geometry']));nv=len(g['vertices']);n=g['physical_spins'];p=saved['edge_exponents'];legs=[[None]*3 for _ in range(nv)]
    for e in g['shared_gauge_edges']:
        for v in e['endpoints']:legs[v['vertex_id']-1][v['local_leg_i']-1]=e['edge_id']-1
    h,basis,_=hamiltonian(legs,saved['nup']);h=h.toarray();lookup={int(b):j for j,b in enumerate(basis)};T=np.eye(len(basis),dtype=complex)
    weights=[0]*n
    for v,es in enumerate(legs):
        c=p[es[0]]%3;q=(p[es[1]]-c)%3
        for a in range(3):weights[3*v+a]=(c+q*a)%3
        rotation=np.zeros_like(h)
        for col,old in enumerate(basis):
            old=int(old);localin=(old>>(3*v))&7;outside=old&~(7<<(3*v))
            for out in range(8):
                if abs(B[out,localin])>1e-14:
                    rotation[lookup[outside|(out<<(3*v))],col]=B[out,localin]
        T=rotation@T
    for e,pe in enumerate(p):weights[3*nv+e]=pe%3
    rotated=T.conj().T@h@T;construct=np.zeros_like(h)
    for v,es in enumerate(legs):
        for i,e in enumerate(es):
            gauge=3*nv+e
            for term in local['rotated_exchange_operators'][i]['nonzero_matrix_units']:
                coefficient=complex(term['real'],term['imag']);localin=term['input'];localout=term['out']
                for col,old in enumerate(basis):
                    old=int(old)
                    if (old>>gauge)&1 or ((old>>(3*v))&7)!=localin:continue
                    new=(old&~(7<<(3*v)))|(localout<<(3*v))|(1<<gauge);row=lookup[new]
                    construct[row,col]-=coefficient;construct[col,row]-=coefficient.conjugate()
    error=np.linalg.norm(rotated-construct);unitarity=np.linalg.norm(T.conj().T@T-np.eye(len(basis)))
    charges=np.array([sum(w*((int(b)>>j)&1) for j,w in enumerate(weights))%3 for b in basis])
    mixing=np.linalg.norm(construct*(charges[:,None]!=charges[None,:]));assert max(error,unitarity,mixing)<1e-11
    sectors=[]
    for charge in range(3):
        chosen=np.flatnonzero(charges==charge);vals=np.linalg.eigvalsh(construct[np.ix_(chosen,chosen)])
        target=next(s['energy'] for s in saved['sectors'] if s['charge']==charge)
        assert abs(vals[0]-target)<1e-11
        sectors.append(dict(charge=charge,dimension=len(chosen),first_eigenvalues=vals[:10].tolist(),
            ground_multiplicity=int(sum(abs(vals-vals[0])<1e-9)),energy_error_vs_original_ed=float(vals[0]-target)))
    records.append(dict(family=saved['family'],geometry=saved['geometry'],nup=saved['nup'],physical_spins=n,
        matter_states_per_star=8,site_charge_weights_physical_order=weights,number_basis_dimension=len(basis),
        exact_full_block_basis_equivalence_error=float(error),unitarity_error=float(unitarity),sector_mixing_error=float(mixing),sectors=sectors))
    print(saved['family'],'exact rotated Hamiltonian and all three onsite charge blocks passed',flush=True)
r=dict(records=records,runtime_seconds=time.time()-start,
    audit=provenance(None,'Independent microscopic ED matrix basis conjugation',dict(tolerance=1e-11,local_roundoff_cutoff=1e-14),'Complete previously validated N_up=2 blocks',['N_up','microscopic winding charge']),
    interpretation='All microscopic states and shared physical gauge spins retained. Direct onsite winding QNs are an exact change of basis, not an effective gauge/clock approximation. Production library implementation still needs validation.')
pathlib.Path('results/charge_basis_ed_audit.json').write_text(json.dumps(r,indent=2)+'\n')
