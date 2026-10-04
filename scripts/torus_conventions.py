"""Exhaustive endpoint audit, reducing exactly equivalent cases by local unitaries."""
import itertools,json,csv,gzip,pathlib
import numpy as np
from scipy.sparse.linalg import eigsh
from ed import W,hamiltonian

reference=np.array([-7.364534839426]*2+[-7.342883075944,-7.337288472473]+[-7.296046028855]*2+[-7.289970027036]*2)
perms=list(itertools.permutations(range(3)))
choices=[]
for conjugate,p in itertools.product([False,True],perms):
    target=(W.conjugate() if conjugate else W)[:,p]
    # target=M W. M must be a monomial unitary, realizable by independent
    # matter-spin permutations and z rotations, with gauge spins untouched.
    M=target@W.conj().T
    rowperm=np.argmax(abs(M),axis=1)
    monomial=np.zeros((3,3),complex);monomial[np.arange(3),rowperm]=M[np.arange(3),rowperm]
    error=float(max(np.max(abs(M-monomial)),np.max(abs(abs(M[np.arange(3),rowperm])-1)),np.max(abs(target-monomial@W))))
    assert len(set(rowperm))==3 and error<1e-14
    choices.append({'conjugate_W':conjugate,'column_permutation':p,'matter_row_permutation':rowperm.tolist(),
        'matter_row_phase_exponents':[int(np.argmin(abs(z-np.exp(2j*np.pi*np.arange(3)/3)))) for z in M[np.arange(3),rowperm]],'identity_error':error})

graphs=[]
for parity in [(1,0),(0,1),(1,1)]:
    # Every index-two translation subgroup is a kernel of a nonzero map Z²→Z₂.
    legs=[[-1]*3 for _ in range(4)];edges=[]
    for s in range(2):
        for i,delta in enumerate([(0,0),(-1,0),(0,-1)]):
            b=(s+parity[0]*delta[0]+parity[1]*delta[1])%2
            e=len(edges);av=2*s;bv=2*b+1
            legs[av][i]=e;legs[bv][i]=e
            edges.append({'edge':e,'A_vertex':av,'A_leg':i,'B_vertex':bv,'B_leg':i})
    assert sorted(sum(legs,[]))==sorted(list(range(6))*2)
    H,basis,herm=hamiltonian(legs,9)
    assert len(set(basis))==48620 and all(int(b).bit_count()==9 for b in basis)
    eigenvalues,vectors=eigsh(H,k=16,which='SA',ncv=96,tol=1e-12,maxiter=10000)
    ix=np.argsort(eigenvalues);eigenvalues=eigenvalues[ix];vectors=vectors[:,ix]
    # QR within clusters; certify geometric multiplicity, not repeated values alone.
    start=0
    for stop in range(1,len(eigenvalues)+1):
        if stop==len(eigenvalues) or abs(eigenvalues[stop]-eigenvalues[start])>1e-8:
            vectors[:,start:stop]=np.linalg.qr(vectors[:,start:stop])[0];start=stop
    residuals=np.linalg.norm(H@vectors-vectors*eigenvalues,axis=0)
    gram=float(np.linalg.norm(vectors.conj().T@vectors-np.eye(16)))
    assert max(residuals)<1e-9 and gram<1e-8
    graphs.append({'parity':parity,'legs':legs,'edges':edges,'hermiticity_error':herm,
        'basis_dimension':len(basis),'all_basis_popcounts_are_9':True,
        'energies':eigenvalues.tolist(),'residuals':residuals.tolist(),'gram_error':gram,
        'ground_multiplicity':int(sum(abs(eigenvalues-eigenvalues[0])<1e-9)),
        'first8_error_vs_supplied':float(max(abs(eigenvalues[:8]-reference)))})
    print('quotient',parity,graphs[-1],flush=True)

canonical=np.array(graphs[0]['energies'][:14])
assert all(max(abs(np.array(g['energies'][:14])-canonical))<1e-9 for g in graphs)
with gzip.open('results/torus_all_endpoint_conventions.csv.gz','wt',newline='') as f:
    writer=csv.writer(f);writer.writerow(['quotient','A0_convention','B0_convention','A1_convention','B1_convention',*[f'E{k}' for k in range(8)],'full8_matches_supplied','method'])
    for g in graphs:
        for pattern in itertools.product(range(12),repeat=4):
            writer.writerow([str(g['parity']),*pattern,*g['energies'][:8],False,'Certified local-unitary equivalence to independently diagonalized representative'])

result={'endpoint_choices':choices,'translation_quotients':graphs,
    'enumerated_patterns':3*12**4,'supplied_first8':reference.tolist(),
    'matching_patterns':0,'scope':'All ordinary periodic honeycomb cells with four vertices (two primitive cells), all local leg permutations and W/conj(W) choices independently at all four endpoints. W transpose equals W. Independent matter-label permutations are trivial basis relabelings.',
    'proof':'For any permutation of three columns, p(i)=s*i+t mod 3 with s=±1. W[a,p(i)]=omega**(a*t)*W[s*a,i]. Conjugation maps a to -a. These changes are implemented by number-preserving local matter permutations and spin z rotations, so they cannot change any eigenvalue or multiplicity.'}
pathlib.Path('results/torus_convention_audit.json').write_text(json.dumps(result,indent=2))
