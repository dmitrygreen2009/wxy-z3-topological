"""Test microscopic string proposals on cached torus states; no eigensolve."""
import json,pathlib,time
import numpy as np
from scipy.sparse import coo_matrix
from ed import hamiltonian,W
start=time.time();omega=np.exp(2j*np.pi/3);legs=[[0,1,2],[0,4,2],[3,4,5],[3,1,5]]
H,basis,_=hamiltonian(legs,9);lookup={int(b):i for i,b in enumerate(basis)};cache=np.load('results/physics48_torus_joint_vectors.npz');G=cache['vectors'];audit=json.load(open('results/physics48_torus.json'))
assert np.array_equal(cache['basis'],basis)
def unitary(p):
 mapped=basis.copy();exponent=np.zeros(len(basis),dtype=np.int64)
 for v,es in enumerate(legs):
  c=p[es[0]];q=(p[es[1]]-c)%3;mapped &= ~(7<<(3*v))
  for a in range(3):
   bit=(basis>>(3*v+a))&1;mapped |= bit<<(3*v+(a+q)%3);exponent+=c*bit
 for e,pe in enumerate(p):exponent+=pe*((basis>>(12+e))&1)
 return coo_matrix((omega**(exponent%3),([lookup[int(x)] for x in mapped],np.arange(len(basis)))),shape=H.shape).tocsr()
def component(vertices):
 rows=[];cols=[];values=[]
 for v in vertices:
  for a in range(3):
   mb=(basis>>(3*v+a))&1
   for i,e in enumerate(legs[v]):
    eb=(basis>>(12+e))&1;ix=np.where(mb!=eb)[0];dst=basis[ix]^(1<<(3*v+a))^(1<<(12+e));rows.extend(lookup[int(x)] for x in dst);cols.extend(ix);values.extend(np.where(mb[ix],-W[a,i],-W[a,i].conjugate()))
 return coo_matrix((values,(rows,cols)),shape=H.shape).tocsr()
def err(M):return float(max(abs(M.data),default=0))
cycles={k:unitary(p) for k,p in audit['edge_exponents'].items()};records=[]
for edges,p in [([0],[1,0,0,0,0,0]),([0,4],[1,0,0,0,2,0]),([0,4,3],[1,0,0,1,2,0])]:
 U=unitary(p);ends=[v for v,row in enumerate(legs) if sum(p[e] for e in row)%3];interior=component([v for v in range(4) if v not in ends]);error=err(U.conj().T@interior@U-interior);assert error<1e-12
 comm={k:err(U@V-V@U) for k,V in cycles.items()};assert max(comm.values())<1e-12
 energies=[]
 for k in range(4):
  state=U@G[:,k];energy=float(np.vdot(state,H@state).real);energies.append(energy-audit['energy'])
 records.append(dict(path_edges_zero_based=edges,edge_exponents=p,endpoint_vertices_zero_based=ends,interior_H_invariance_error=error,CGS_commutator_errors=comm,energy_costs_four_ground_states=energies,endpoint_operator_norm_bound=2*np.sqrt(3)*len(ends),interpretation='Exact unitary with matter action; H changes only at path endpoints, but all CGS cycle charges are unchanged. Bounded energy is automatic from the local symmetry, not a deconfinement test. Path length is not shortest endpoint separation on this tiny quotient.'))
# A physically number-neutral ladder pair changes *symmetry* charges.
# Matter permutations in CGS commute with these gauge ladders, so no missing
# matter dressing is required for this covariance identity.
raise_bit=12;lower_bit=15;ix=[i for i,b in enumerate(basis) if not (int(b)>>raise_bit&1) and (int(b)>>lower_bit&1)]
O=coo_matrix((np.ones(len(ix)),([lookup[int(basis[i])^(1<<raise_bit)^(1<<lower_bit)] for i in ix],ix)),shape=H.shape).tocsr()
shifts={name:(p[0]-p[3])%3 for name,p in audit['edge_exponents'].items()};cov={name:err(U@O-omega**shifts[name]*O@U) for name,U in cycles.items()};assert max(cov.values())<1e-12
k=next(k for k,r in enumerate(audit['joint_cycle_states']) if (r['Y0_q'],r['Y1_q'])==(1,1));state=O@G[:,k];weight=float(np.vdot(state,state).real);assert weight>1e-12;state/=np.sqrt(weight);energy=float(np.vdot(state,H@state).real)
charge={name:[float(np.vdot(state,U@state).real),float(np.vdot(state,U@state).imag)] for name,U in cycles.items()}
r=dict(open_CGS_unitaries=records,neutral_ladder_pair=dict(operator='sigma^+_edge0 sigma^-_edge3',number_change=0,covariance_shifts=shifts,covariance_errors=cov,ground_initial_charges=[1,1],final_charges=[2,0],state_norm_before_normalization=weight,energy=energy,energy_cost=energy-audit['energy'],energy_residual=float(np.linalg.norm(H@state-energy*state)),exact_sector_minimum_reused=-7.296046028854623,exact_sector_minimum_cost=-7.296046028854623-audit['energy'],cycle_expectations=charge,interpretation='CGS-charge-covariant, not invariant under each local/cycle symmetry. Creates two winding-charge changes on this quotient; not a demonstrated localized anyon pair. No path or separation dependence can be inferred.'),
 definition_for_future_test='For independently defined contractible CGS generators B_f, project onto a fixed reference q_f pattern except two endpoints q_f shifted +1,-1. E_pair(R) is the minimum of unchanged H in the nonzero joint projector range at specified N_up, boundary and winding sectors; subtract the compatible defect-free minimum. This tests interaction of these symmetry defects only after their identification as the desired physical excitation. If endpoints can move within a background with varying q_f, reference pattern must first be resolved.',
 no_confinement_fit=True,reason='Tiny quotient lacks independent distant plaquettes. Saved larger MPS are unconverged and not separately optimized in two-defect sectors. Applying a trial string gives only an upper bound, not lowest defect energy. Closed CGS loop expectations and matter-dressed open symmetry strings cannot alone diagnose deconfinement.',runtime_seconds=time.time()-start,solver='No eigensolver: sparse operator identities and cached-vector expectations')
r['analysis_git_commit']=__import__('subprocess').check_output(['git','rev-parse','HEAD'],text=True).strip()
pathlib.Path('results/physics48_strings.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
