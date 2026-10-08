"""New small symmetry-block excitations only; reuse and deflate saved ground vectors."""
import pathlib,json,itertools,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh,LinearOperator
from ed import hamiltonian
start=time.time();omega=np.exp(2j*np.pi/3);legs=[[0,1,2],[0,4,2],[3,4,5],[3,1,5]]
H,basis,_=hamiltonian(legs,9);lookup={int(b):i for i,b in enumerate(basis)}
cache=np.load('results/physics48_torus_joint_vectors.npz');G=cache['vectors'];audit=json.load(open('results/physics48_torus.json'))
# Recover the exact two commuting monomial generators.
maps=[];phases=[]
for p in [[1,0,2,0,0,0],[0,0,0,1,0,2]]:
 dst=[];phase=[]
 for old in basis:
  old=int(old);new=old;exp=0
  for v,es in enumerate(legs):
   c=p[es[0]];q=(p[es[1]]-c)%3;new&=~(7<<(3*v))
   for a in range(3):
    bit=(old>>(3*v+a))&1;new|=bit<<(3*v+(a+q)%3);exp+=c*bit
  exp+=sum(p[e]*((old>>(12+e))&1) for e in range(6))
  dst.append(lookup[new]);phase.append(omega**(exp%3))
 maps.append(np.array(dst));phases.append(np.array(phase))
seen=np.zeros(len(basis),bool);orbits=[]
for rep in range(len(basis)):
 if seen[rep]:continue
 terms=[]
 for a,b in itertools.product(range(3),repeat=2):
  v=rep;z=1+0j
  for k,times in [(0,a),(1,b)]:
   for _ in range(times):z*=phases[k][v];v=maps[k][v]
  terms.append((a,b,v,z));seen[v]=True
 orbits.append(terms)
records=[]
for charge in [(0,0),(1,0),(1,1)]:
 rows=[];cols=[];data=[];dim=0
 for terms in orbits:
  amps={}
  for a,b,v,z in terms:amps[v]=amps.get(v,0)+z*omega**((-charge[0]*a-charge[1]*b)%3)
  norm=np.sqrt(sum(abs(z)**2 for z in amps.values()))
  if norm<1e-10:continue
  for v,z in amps.items():
   if abs(z)>1e-12:rows.append(v);cols.append(dim);data.append(z/norm)
  dim+=1
 P=coo_matrix((data,(rows,cols)),shape=(len(basis),dim)).tocsr();hb=(P.conj().T@H@P).tocsr()
 assert np.max(abs((P.conj().T@P).diagonal()-1))<1e-12
 ix=[k for k,r in enumerate(audit['joint_cycle_states']) if (r['Y0_q'],r['Y1_q'])==charge]
 shift=20.;g=P.conj().T@G[:,ix] if ix else np.zeros((dim,0),complex)
 if ix:
  assert np.linalg.norm(g.conj().T@g-np.eye(len(ix)))<1e-10
  assert np.linalg.norm(hb@g-audit['energy']*g)<1e-9
 op=LinearOperator(hb.shape,matvec=lambda x:hb@x+shift*g@(g.conj().T@x),dtype=complex)
 rng=np.random.default_rng(7310+charge[0]*3+charge[1]);x=rng.normal(size=dim)+1j*rng.normal(size=dim)
 if ix:x-=g@(g.conj().T@x)
 t=time.time();e,V=eigsh(op,k=2,which='SA',tol=1e-12,ncv=60,maxiter=10000,v0=x);order=np.argsort(e);e=e[order];V=V[:,order]
 residual=np.linalg.norm(hb@V-V*e,axis=0);assert max(residual)<1e-9
 # Independent full-H residual and charge eigenvalues on recovered eigenvectors.
 full=P@V;assert max(np.linalg.norm(H@full-full*e,axis=0))<1e-9
 records.append(dict(charges=charge,dimension=dim,energies=e.tolist(),residuals=residual.tolist(),saved_ground_vectors_deflated=len(ix),runtime_seconds=time.time()-t))
 print(records[-1],flush=True)
 pathlib.Path('results/physics48_torus_excitations.json').write_text(json.dumps(dict(records=records,seed=7310,tolerance=1e-12,
  method='Exact CGS orbit blocks; cached ground eigenvectors deflated, no ground-state Lanczos repeated',excited_vectors_previously_available=False,
  symmetry_equivalences='R0 and R1 invert separate charges; Tx exchanges charges. Three charge orbits cover all nine sectors.',runtime_seconds=time.time()-start),indent=2)+'\n')
 np.savez_compressed(f'results/physics48_torus_excited_q{charge[0]}{charge[1]}.npz',basis=basis,energies=e,vectors=full)
