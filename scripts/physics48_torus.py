"""Read-only saved-groundspace symmetry extraction; no eigensolver on H."""
import itertools,json,pathlib,time,sys
import numpy as np
from scipy.sparse import coo_matrix
from ed import hamiltonian
started=time.time();omega=np.exp(2j*np.pi/3)
legs=[[0,1,2],[0,4,2],[3,4,5],[3,1,5]]
H,basis,herm=hamiltonian(legs,9)
d=np.load('results/torus_groundspace.npz');V=d['ground_vectors']
assert np.array_equal(basis,d['basis']);assert V.shape==(48620,4)
E=-7.364534839425767;res=np.linalg.norm(H@V-E*V,axis=0);gram=np.linalg.norm(V.conj().T@V-np.eye(4))
assert max(res)<1e-9 and gram<1e-8
lookup={int(b):i for i,b in enumerate(basis)}

def monomial(mapping,phases=None):
 dst=np.array([lookup[mapping(int(b))] for b in basis]);phase=np.ones(len(basis),complex) if phases is None else phases
 return coo_matrix((phase,(dst,np.arange(len(basis)))),shape=H.shape).tocsr()

def cycle(p):
 assert all(sum(p[e] for e in row)%3==0 for row in legs)
 dst=[];phase=[]
 for b in basis:
  old=int(b);new=old;exponent=0
  for v,edges in enumerate(legs):
   c=p[edges[0]];q=(p[edges[1]]-c)%3;new &= ~(7<<(3*v))
   for a in range(3):
    bit=(old>>(3*v+a))&1;new|=bit<<(3*v+(a+q)%3);exponent+=c*bit
  exponent+=sum(p[e]*((old>>(12+e))&1) for e in range(6))
  dst.append(lookup[new]);phase.append(omega**(exponent%3))
 return coo_matrix((phase,(dst,np.arange(len(basis)))),shape=H.shape).tocsr()

def perm(bits):
 return lambda b:sum(((b>>i)&1)<<bits[i] for i in range(18))
ops={};ps={'Y0':[1,0,2,0,0,0],'Y1':[0,0,0,1,0,2],'X':[1,2,0,1,2,0],'P_hex':[1,0,2,2,0,1]}
for name,p in ps.items():ops[name]=cycle(p)
# Full index-two primitive translation x, vertices A0,B0,A1,B1.
tbits=list(range(6,12))+list(range(0,6))+list(range(15,18))+list(range(12,15))
ops['Tx']=monomial(perm(tbits))
# Exact unitary charge conjugation: full spin flip and exchange matter rows1,2.
pbits=list(range(18))
for v in range(4):pbits[3*v+1],pbits[3*v+2]=pbits[3*v+2],pbits[3*v+1]
ops['F']=monomial(lambda b:perm(pbits)(((1<<18)-1)^b))
# Independent reflections of each identified parallel-edge pair, with the
# exact local matter phases/permutations implied by W column relabeling.
for cell in range(2):
 bits=list(range(18))
 for v in [2*cell,2*cell+1]:bits[3*v+1],bits[3*v+2]=bits[3*v+2],bits[3*v+1]
 e0,e2=12+3*cell,14+3*cell;bits[e0],bits[e2]=bits[e2],bits[e0]
 phases=np.array([omega**(sum(2*a*((int(b)>>(3*v+a))&1) for v in [2*cell,2*cell+1] for a in range(3))%3) for b in basis])
 ops['R'+str(cell)]=monomial(perm(bits),phases)
records={};sub={}
for name,U in ops.items():
 comm=H@U-U@H;error=max(abs(comm.data),default=0);assert error<1e-13,(name,error)
 G=V.conj().T@(U@V);leak=np.linalg.norm(U@V-V@G);assert leak<1e-8
 sub[name]=G
 records[name]=dict(full_commutator_max=float(error),groundspace_invariance_residual=float(leak),matrix_real=G.real.tolist(),matrix_imag=G.imag.tolist(),eigenvalues=[[float(z.real),float(z.imag)] for z in np.linalg.eigvals(G)])
# Cycles commute; diagonalize a generic Hermitian joint label.
K=sum((i+1)*(U+U.conj().T)/2+((i+1)**2)*((U-U.conj().T)/(2j)) for i,U in enumerate([sub['Y0'],sub['X']]))
_,R=np.linalg.eigh(K);G=V@R
labels=[]
for k in range(4):
 row={'state':k,'energy':float(np.vdot(G[:,k],H@G[:,k]).real),'charges':{},'unitary_variances':{}}
 for name,U in ops.items():
  z=np.vdot(G[:,k],U@G[:,k]);row['charges'][name]=[float(z.real),float(z.imag)];row['unitary_variances'][name]=float(max(0,1-abs(z)**2))
  if name in ps:row[name+'_q']=int(np.argmin(abs(z-omega**np.arange(3))))
 labels.append(row)
reflection_connections={}
for name in ['R0','R1','Tx','F']:
 M=G.conj().T@(ops[name]@G);reflection_connections[name]=np.abs(M).tolist()
I=coo_matrix((np.ones(len(basis)),(np.arange(len(basis)),np.arange(len(basis)))),shape=H.shape).tocsr()
def maxentry(M):return float(max(abs(M.data),default=0))
relations={
 'R0_squared_identity':maxentry(ops['R0']@ops['R0']-I),
 'R1_squared_identity':maxentry(ops['R1']@ops['R1']-I),
 'R0_R1_commute':maxentry(ops['R0']@ops['R1']-ops['R1']@ops['R0']),
 'R0_inverts_Y0':maxentry(ops['R0']@ops['Y0']@ops['R0']-ops['Y0'].getH()),
 'R1_inverts_Y1':maxentry(ops['R1']@ops['Y1']@ops['R1']-ops['Y1'].getH()),
 'R0_commutes_Y1':maxentry(ops['R0']@ops['Y1']-ops['Y1']@ops['R0']),
 'R1_commutes_Y0':maxentry(ops['R1']@ops['Y0']-ops['Y0']@ops['R1']),
 'Y0_Y1_X_identity':maxentry(ops['Y0']@ops['Y1']@ops['X']-I),
 'P_hex_equals_Y0_Y1adjoint':maxentry(ops['P_hex']-ops['Y0']@ops['Y1'].getH())}
assert max(relations.values())<1e-13
commuting={f'{a},{b}':float(np.linalg.norm(sub[a]@sub[b]-sub[b]@sub[a])) for a,b in itertools.combinations(ops,2)}
out=dict(energy=E,ground_vector_residuals=res.tolist(),gram_error=float(gram),operators=records,joint_cycle_states=labels,restricted_commutators=commuting,edge_exponents=ps,groundspace_symmetry_connections=reflection_connections,exact_group_relations=relations,
 exact_groundspace_cycle_relation='Y0*Y1*X = identity because exponent sum is 2 on all edges and N_up=9',
 runtime_seconds=time.time()-started,diagonalizations_of_H_performed=False)
pathlib.Path('results/physics48_torus.json').write_text(json.dumps(out,indent=2)+'\n')
np.savez_compressed('results/physics48_torus_joint_vectors.npz',basis=basis,vectors=G)
print(json.dumps({'labels':labels,'commutators':commuting,'residuals':res.tolist()},indent=2),flush=True)
