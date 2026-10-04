"""Exact extra microscopic reflection symmetry of the narrow zigzag quotient.
Validate the operator identities on every number block of the full ten-spin
parallel-edge pair; no clock or gauge approximation is made.
"""
import json,pathlib,time
import numpy as np
from scipy.sparse import coo_matrix,eye
from ed import hamiltonian
from provenance import provenance
omega=np.exp(2j*np.pi/3)
geometry_path="geometry/zigzag_L1_w1_star.json"
geometry=json.loads(pathlib.Path(geometry_path).read_text())
legs=[[int(e)-1 for e in row] for row in geometry["endpoint_leg_table"]]
assert geometry["physical_spins"]==10 and legs==[[0,1,2],[0,3,2]]
def operators(basis):
    lookup={int(b):k for k,b in enumerate(basis)}
    reflection=basis.copy();phase=np.zeros(len(basis),dtype=int)
    cycle=basis.copy();cyclephase=np.zeros(len(basis),dtype=int)
    for v in range(2):
        reflection &= ~(7<<(3*v));cycle &= ~(7<<(3*v))
        for a in range(3):
            bit=(basis>>(3*v+a))&1
            reflection |= bit<<(3*v+(-a)%3);phase-=a*bit
            cycle |= bit<<(3*v+(a+2)%3);cyclephase+=bit
    reflection &= ~((1<<6)|(1<<8))
    reflection |= ((basis>>6)&1)<<8;reflection |= ((basis>>8)&1)<<6
    cyclephase+=((basis>>6)&1)+2*((basis>>8)&1)
    def matrix(mapped,exponents):
        return coo_matrix((omega**(exponents%3),([lookup[int(b)] for b in mapped],np.arange(len(basis)))),shape=(len(basis),len(basis))).tocsr()
    return matrix(reflection,phase),matrix(cycle,cyclephase)
def error(A):return float(max(abs(A.data),default=0))
points=[]
for nup in range(11):
    started=time.time();H,basis,herm=hamiltonian(legs,nup);R,U=operators(basis);I=eye(len(basis))
    errors=dict(reflection_commutator=error(H@R-R@H),cycle_commutator=error(H@U-U@H),
        reflection_square=error(R@R-I),cycle_cube=error(U@U@U-I),dihedral_relation=error(R@U@R-U.getH()))
    assert max(errors.values())<1e-13
    energies,vectors=np.linalg.eigh(H.toarray());ground=vectors[:,abs(energies-energies[0])<1e-9]
    charge_eigenvalues=np.linalg.eigvals(ground.conj().T@(U@ground))
    charges=[int(np.argmin(abs(z-omega**np.arange(3)))) for z in charge_eigenvalues]
    assert max(abs(z-omega**q) for z,q in zip(charge_eigenvalues,charges))<1e-9
    points.append(dict(nup=nup,dimension=len(basis),ground_energy=float(energies[0]),ground_multiplicity=ground.shape[1],
        ground_cycle_charges=charges,operator_errors=errors,runtime_seconds=time.time()-started))
output=dict(geometry_file=geometry_path,legs_zero_based=legs,matter_spins=[list(range(3)),list(range(3,6))],gauge_spins=list(range(6,10)),
    parallel_physical_edges=[0,2],reflection='Exchange physical gauge edges 0 and 2; at each endpoint map matter species a to -a mod 3 and multiply each ket by omega^(-sum_a a*n_up[a]).',
    points=points,audit=provenance(None,'Complete-number-block exact operator identities',{'identity_tolerance':1e-13},'All 1024 physical kets across eleven exact number blocks',['N_up']),
    interpretation='Parallel edges create an additional local S3 symmetry: R^2=U^3=I and RUR=U†. Nontrivial cycle-charge doublets are exactly degenerate when present. This quotient-specific enhancement is not a proof of a bulk phase or of extensive ground degeneracy.')
pathlib.Path('results/parallel_edge_symmetry_audit.json').write_text(json.dumps(output,indent=2))
print('Passed full-Hilbert local S3 identities in all eleven number blocks.')
