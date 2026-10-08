"""New analytic bulk trial and small-matrix operator bound; no optimization of MPS."""
import json,math,pathlib
import numpy as np
from scipy.optimize import minimize_scalar
from ed import hamiltonian
omega=np.exp(2j*np.pi/3);W=np.array([[1,1,1],[1,omega,omega.conjugate()],[1,omega.conjugate(),omega]])/np.sqrt(3)
g=np.array([1,1,omega]);m=(W@g).conjugate()
assert np.max(abs(abs(m)-1))<1e-14
checks=[]
for name,legs,nup in [('star',[[0,1,2]],3),('pair',[[0,1,2],[0,3,4]],5),('torus',[[0,1,2],[0,4,2],[3,4,5],[3,1,5]],9)]:
 h,basis,herm=hamiltonian(legs,nup);nv=len(legs);ns=3*nv+max(sum(legs,[]))+1
 z=np.empty(ns,complex);z[:3*nv]=np.tile(m,nv);assigned={}
 for edges in legs:
  for i,e in enumerate(edges):
   if e in assigned:assert abs(assigned[e]-g[i].conjugate())<1e-14
   assigned[e]=g[i].conjugate();z[3*nv+e]=g[i].conjugate()
 state=np.array([np.prod([z[j] for j in range(ns) if int(b)>>j&1]) for b in basis]);state/=np.linalg.norm(state)
 energy=float(np.vdot(state,h@state).real);prediction=-6*nv*nup*(ns-nup)/(ns*(ns-1))
 assert abs(energy-prediction)<1e-12
 checks.append(dict(geometry=name,nup=nup,spins=ns,energy=energy,prediction=prediction,error=abs(energy-prediction),hamiltonian_hermiticity=herm))
# Exact local operator inequality h_v >= epsilon(lambda)+lambda*(Nm+Ng/2).
# Shared gauge number then sums exactly once over stars, including matter once.
blocks=[]
for n in range(7):
 h,b,_=hamiltonian([[0,1,2]],n);q=np.array([(int(x)&7).bit_count()+.5*(int(x)>>3).bit_count() for x in b]);blocks.append((h.toarray(),q))
def epsilon(lam):return min(float(np.linalg.eigvalsh(h-lam*np.diag(q))[0]) for h,q in blocks)
def threshold(lam):return (-3-2*(epsilon(lam)-1e-10))/(9*lam)
opt=minimize_scalar(lambda x:-threshold(x),bounds=(-4,-.2),method='bounded',options={'xatol':1e-10});lam=float(opt.x)
# Certify full eigen-decompositions at the single chosen multiplier, using
# Frobenius reconstruction error and orthogonality error; conservative margin.
certs=[];lower=[]
for h,q in blocks:
 A=h-lam*np.diag(q);evals,V=np.linalg.eigh(A);rec=np.linalg.norm(A-(V*evals)@V.conj().T);orth=np.linalg.norm(V.conj().T@V-np.eye(len(q)))
 error=rec+np.linalg.norm(A)*orth+1e-10
 lower.append(float(evals[0]-error));certs.append(dict(dimension=len(q),reconstruction_error=float(rec),orthogonality_error=float(orth),conservative_margin=float(error)))
eps=min(lower);rho=(-3-2*eps)/(9*lam)
r=dict(trial_gauge_plus_phases=[[z.real,z.imag] for z in g],trial_matter_ket_up_phases=[[z.real,z.imag] for z in m],trial_energy_per_primitive_cell=-3,trial_mean_filling=.5,trial_breaks_U1=True,trial_is_valid_without_assuming_ground_filling=True,fixed_number_benchmark_expectations=checks,
 analytic_density_window=[1/(3*math.sqrt(2)),1-1/(3*math.sqrt(2))],
 cluster_inequality='H>=Nv*epsilon(lambda)+lambda*N_up, epsilon=min eig(h_star-lambda*(Nm+Ng/2)); valid on closed degree-three bulk',cluster_multiplier=lam,cluster_epsilon_conservative=eps,cluster_density_window=[rho,1-rho],small_matrix_certificates=certs,
 limitations='Bulk degree-three lattice only. No narrow-cylinder energy density imported as a 2D trial. Number-projected trial retains the same limiting energy and half density, but the global bound does not assume the ground filling. Small-matrix inequality uses an explicit conservative 1e-10 margin; numerical certificate, not symbolic optimality proof.')
r['particle_hole_locality_theorem']={'statement': 'The thermodynamic canonical ground energy density e(rho) has a minimum at rho=1/2, without assuming an unbroken particle-hole symmetry in a pure phase.', 'proof': 'For any sequence of canonical states at density rho, construct two adjacent macroscopic regions in that state and its exact particle-hole image. Their densities average to 1/2 and their energy densities coincide. Finite-range bounded interactions add only O(boundary area) energy, which vanishes per volume. Thus e(1/2)<=e(rho). Equivalently canonical e is convex by domain mixing and particle-hole symmetric.', 'assumptions': ['thermodynamic energy density exists', 'short-range bounded microscopic interactions', 'unrestricted thermodynamic limit allowing macroscopic domains'], 'not_proven': ['uniqueness of the minimizing density', 'half filling of every pure clustering ground state', 'unbroken particle-hole or primitive translations', 'spectral or charge gap'], 'alternative': 'If distinct rho and 1-rho are global pure-phase minima, convexity gives a flat energy-density interval containing half filling; the half-filled minimizer may be phase coexistence rather than a homogeneous phase.'}
r['analysis_git_commit']=__import__('subprocess').check_output(['git','rev-parse','HEAD'],text=True).strip()
pathlib.Path('results/physics48_density_bound.json').write_text(json.dumps(r,indent=2)+'\n');print('Density windows: analytic',r['analytic_density_window'],'cluster',r['cluster_density_window'],'lambda',lam)
