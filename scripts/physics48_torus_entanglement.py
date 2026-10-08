"""Spatial Schmidt analysis of saved torus states; no new energy optimization."""
import json,pathlib
import numpy as np
x=np.load('results/physics48_torus_joint_vectors.npz');b=x['basis'];V=x['vectors'];lookup={int(v):i for i,v in enumerate(b)}
a=json.load(open('results/physics48_torus.json'));k=next(i for i,r in enumerate(a['joint_cycle_states']) if r['Y0_q']==r['Y1_q']==1)
omega=np.exp(2j*np.pi/3)
def reflect(v,cell):
 bits=list(range(18))
 for star in [2*cell,2*cell+1]:bits[3*star+1],bits[3*star+2]=bits[3*star+2],bits[3*star+1]
 i,j=12+3*cell,14+3*cell;bits[i],bits[j]=bits[j],bits[i]
 y=np.empty_like(v)
 for col,old in enumerate(b):
  old=int(old);new=sum(((old>>i)&1)<<bits[i] for i in range(18))
  exponent=sum(2*species*((old>>(3*star+species))&1) for star in [2*cell,2*cell+1] for species in range(3))%3
  y[lookup[new]]=omega**exponent*v[col]
 return y
v=V[:,k];v/=np.linalg.norm(v);orbit=np.column_stack([v,reflect(v,1),reflect(v,0),reflect(reflect(v,1),0)])
assert np.linalg.norm(orbit.conj().T@orbit-np.eye(4))<1e-10
left=[0,1,2,3,4,5,12,14,16];right=[i for i in range(18) if i not in left]
li=sum(((b>>s)&1)<<i for i,s in enumerate(left));ri=sum(((b>>s)&1)<<i for i,s in enumerate(right))
blocks=[(np.array([i for i in range(512) if i.bit_count()==q]),np.array([i for i in range(512) if i.bit_count()==9-q])) for q in range(10)]
def entropy(v):
 M=np.zeros((512,512),complex);M[li,ri]=v
 p=np.concatenate([np.linalg.svd(M[np.ix_(l,r)],compute_uv=False)**2 for l,r in blocks]);p=p[p>1e-15];p/=sum(p)
 return float(-sum(p*np.log(p))),p
s0,ps=entropy(orbit[:,0]);singles=[entropy(orbit[:,i])[0] for i in range(4)]
assert max(abs(np.array(singles)-s0))<1e-10
rng=np.random.default_rng(7311);checks=[]
for c in [np.array([1,0,0,1])/np.sqrt(2),np.ones(4)/2]+[rng.normal(size=4)+1j*rng.normal(size=4) for _ in range(4)]:
 c=c/np.linalg.norm(c);p=np.linalg.svd(c.reshape(2,2),compute_uv=False)**2;p=p[p>1e-15]
 expected=s0-float(sum(p*np.log(p)));measured,_=entropy(orbit@c);assert abs(expected-measured)<1e-9
 checks.append(dict(actual=measured,predicted=expected,error=abs(expected-measured)))
out=dict(left_physical_sites_zero_based=left,right_physical_sites_zero_based=right,charge_basis_entropies=singles,base_entropy=s0,
 logical_bell_entropy=checks[0]['actual'],entropy_range=[s0,s0+np.log(2)],superposition_checks=checks,baseline_schmidt_probabilities=ps.tolist(),
 explanation='Y0 and R0 have support entirely on the left; Y1 and R1 entirely on the right. The 4D ground space is a left logical doublet times a right logical doublet. S(c)=S_base+S(singular values squared of 2x2 coefficient matrix). Extra ln2 is ordinary logical superposition entanglement, not thermodynamic TEE.',
 original_cached_entropy='results/torus_entropy.json',new_hamiltonian_diagonalizations=False)
pathlib.Path('results/physics48_torus_entanglement.json').write_text(json.dumps(out,indent=2)+'\n');print('Sbase',s0,'Bell',checks[0]['actual'],'maxerror',max(c['error'] for c in checks))
