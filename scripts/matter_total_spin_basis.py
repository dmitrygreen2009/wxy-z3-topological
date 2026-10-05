"""Exact requested SU(2) multiplicity decomposition; no production rerun."""
import json,pathlib,time
import numpy as np
from ed import W
from provenance import provenance

started=time.time();omega=np.exp(2j*np.pi/3)
Brecord=json.load(open('results/matter_charge_basis_audit.json'))
B=np.array(Brecord['basis_real'])+1j*np.array(Brecord['basis_imag'])
bits=[0,1,6,7,2,3,4,5];phases=np.array([1,1,1,1,1,-1,1,-1])
U=B[:,bits]*phases
labels=[dict(S=S,m=m,cycle_charge=k,cycle_eigenvalue_real=float((omega**k).real),
    cycle_eigenvalue_imag=float((omega**k).imag),matter_nup=int(m+1.5),
    existing_mps_mode_bits=b,phase_relative_to_existing_encoder=int(sign))
    for S,m,k,b,sign in zip([1.5]*4+[.5]*4,[-1.5,-.5,.5,1.5,-.5,.5,-.5,.5],
        [0,0,0,0,1,1,2,2],bits,phases)]
raises=[]
for a in range(3):
    op=np.zeros((8,8),complex)
    for old in range(8):
        if not (old>>a)&1:op[old|(1<<a),old]=1
    raises.append(op)
Sp=sum(raises);Sm=Sp.conj().T
Sz=np.diag([old.bit_count()-1.5 for old in range(8)])
S2=Sz@Sz+(Sp@Sm+Sm@Sp)/2
C=np.zeros((8,8),complex)
for old in range(8):
    new=sum(((old>>a)&1)<<((a+1)%3) for a in range(3));C[new,old]=1
expected_raise=np.zeros((8,8),complex)
for lower,upper,coefficient in [(0,1,np.sqrt(3)),(1,2,2),(2,3,np.sqrt(3)),(4,5,1),(6,7,1)]:
    expected_raise[upper,lower]=coefficient
checks=dict(unitarity=np.linalg.norm(U.conj().T@U-np.eye(8)),
    total_spin=np.linalg.norm(U.conj().T@S2@U-np.diag([r['S']*(r['S']+1) for r in labels])),
    magnetization=np.linalg.norm(U.conj().T@Sz@U-np.diag([r['m'] for r in labels])),
    cyclic_permutation=np.linalg.norm(U.conj().T@C@U-np.diag([omega**r['cycle_charge'] for r in labels])),
    spin_ladder_phase=np.linalg.norm(U.conj().T@Sp@U-expected_raise))
H=np.zeros((64,64),complex);transformed=np.zeros_like(H);mode=np.zeros_like(H)
for i in range(3):
    A=sum(W[a,i]*raises[a].conj().T for a in range(3))
    term=np.kron(raises[i],A);H-=term+term.conj().T
    term=np.kron(raises[i],U.conj().T@A@U);transformed-=term+term.conj().T
    term=np.kron(raises[i],B.conj().T@A@B);mode-=term+term.conj().T
V=np.kron(np.eye(8),U);P=B.conj().T@U;packed=np.kron(np.eye(8),P)
checks.update(hamiltonian_matrix=np.linalg.norm(V.conj().T@H@V-transformed),
    existing_mps_matrix=np.linalg.norm(packed.conj().T@mode@packed-transformed),
    full_spectrum=float(np.max(np.abs(np.linalg.eigvalsh(H)-np.linalg.eigvalsh(transformed)))),
    hermiticity=np.linalg.norm(transformed-transformed.conj().T))
assert all(error<1e-12 for error in checks.values()),checks
# S of a local star is a basis label, not an extra Hamiltonian conservation law.
local_spin_commutator=float(np.linalg.norm(H@np.kron(np.eye(8),S2)-np.kron(np.eye(8),S2)@H))
assert local_spin_commutator>1e-6
for row in labels:
    b=row['existing_mps_mode_bits'];assert b.bit_count()==row['matter_nup']
    assert sum(a*((b>>a)&1) for a in range(3))%3==row['cycle_charge']
geometry=dict(computational_rows='Integers 0..7; bit a=1 means matter spin a+1 up',
    transformation_convention='Columns of U are coupled-spin kets in the computational basis; coordinates transform with U dagger',
    labels=labels,unitary_real=U.real.tolist(),unitary_imag=U.imag.tolist(),
    states_retained=8,winding_charge_rule='Sum_v(c_v*N_v+s_v*k_v)+sum_e p_e*n_e mod 3',
    c_v='p of local leg 1',s_v='p(leg 2)-p(leg 1) mod 3',
    individual_site_encoding='Use the saved three mode bits. Their additive weights c_v+s_v*a encode the complete multiplicity charge exactly.',
    local_total_spin_is_conserved=False)
pathlib.Path('geometry/matter_total_spin_basis.json').write_text(json.dumps(geometry,indent=2)+'\n')
result=dict(errors={k:float(v) for k,v in checks.items()},local_spin_squared_commutator=local_spin_commutator,
    full_star_energies=np.linalg.eigvalsh(transformed).tolist(),runtime_seconds=time.time()-started,
    audit=provenance(None,'Independent full 8-state SU2/C3 and full 64-state star operator equality',
        dict(tolerance=1e-12),'All matter states and all gauge states',['Physical total N_up','Exact CGS cycle charge']),
    production_relation='Exactly the existing full-state QN basis with a signed column permutation. No solver or checkpoint change required.',
    classification='REJECTED / DIAGNOSTIC ONLY',
    interpretation='Exact representation validation; existing validated MPS/ED/projector comparisons are preserved, not rerun. No phase inference from this algebra.')
pathlib.Path('results/matter_total_spin_basis_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print('Exact quartet plus two cycle-labeled doublets: all matrix checks passed.',checks,flush=True)
