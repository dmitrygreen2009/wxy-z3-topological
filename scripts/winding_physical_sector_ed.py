"""New independent ED of fixed microscopic winding blocks at audited fillings.
This resolves sectors not identified by the earlier unrestricted spectra;
it does not overwrite or rerun the torus or original cylinder benchmarks.
The complete eight-state identity is proved in real_charge_hamiltonian_audit.
"""
import json, pathlib, time, sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from provenance import provenance

fixture='--fixture' in sys.argv
started=time.time();records=[];seed=7258;rng=np.random.default_rng(seed)
for family,nup in [('zigzag',2 if fixture else 9),('armchair',2 if fixture else 8)]:
    geometry=f'geometry/{family}_L2_w1_star.json';g=json.load(open(geometry))
    nv=len(g['vertices']);n=g['physical_spins'];legs=[[None]*3 for _ in range(nv)]
    for e in g['shared_gauge_edges']:
        for endpoint in e['endpoints']:
            legs[endpoint['vertex_id']-1][endpoint['local_leg_i']-1]=e['edge_id']-1
    cycle=json.load(open(f'geometry/cgs_cycles/{family}_L2_w1_star.json'))['cycles'][0]
    p=cycle['edge_exponents'];weights=np.zeros(n,dtype=int)
    for v,es in enumerate(legs):
        c=p[es[0]]%3;shift=(p[es[1]]-c)%3
        assert all((c+shift*i)%3==p[e]%3 for i,e in enumerate(es))
        for a in range(3):weights[3*v+a]=(c+shift*a)%3
    weights[3*nv:]=p
    allbits=np.arange(1<<n,dtype=np.int64)
    basis=allbits[np.bitwise_count(allbits.astype(np.uint64))==nup];del allbits
    charges=sum(int(w)*((basis>>j)&1) for j,w in enumerate(weights))%3
    sectors=[]
    for q in range(3):
        t=time.time();selected=basis[charges==q];rows=[];cols=[];values=[]
        for v,es in enumerate(legs):
            for i,e in enumerate(es):
                m=3*v+i;gauge=3*nv+e;j,k=[3*v+a for a in range(3) if a!=i]
                valid=(((selected>>m)&1)==1)&(((selected>>gauge)&1)==0)
                col=np.flatnonzero(valid);old=selected[col];new=old^((1<<m)|(1<<gauge))
                row=np.searchsorted(selected,new);assert np.array_equal(selected[row],new)
                same=(((old>>j)&1)==((old>>k)&1))
                val=np.where(same,-1.,1/np.sqrt(3))
                rows.append(row);cols.append(col);values.append(val)
                valid=(((selected>>m)&1)==0)&(((selected>>j)&1)==1)&(((selected>>k)&1)==1)&(((selected>>gauge)&1)==0)
                col=np.flatnonzero(valid);old=selected[col]
                new=old^((1<<m)|(1<<j)|(1<<k)|(1<<gauge))
                row=np.searchsorted(selected,new);assert np.array_equal(selected[row],new)
                rows.append(row);cols.append(col);values.append(np.full(len(col),-2/np.sqrt(3)))
        directed=coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(len(selected),len(selected))).tocsr()
        H=directed+directed.T;del directed,rows,cols,values
        assert (H-H.T).nnz==0
        evals,vecs=eigsh(H,k=8,which='SA',ncv=80,tol=1e-12,maxiter=10000,v0=rng.normal(size=len(selected)))
        ix=np.argsort(evals);evals=evals[ix];vecs=vecs[:,ix]
        first=0
        for stop in range(1,len(evals)+1):
            if stop==len(evals) or abs(evals[stop]-evals[first])>1e-8:
                vecs[:,first:stop]=np.linalg.qr(vecs[:,first:stop])[0];first=stop
        residuals=np.linalg.norm(H@vecs-vecs*evals,axis=0)
        gram=float(np.linalg.norm(vecs.T@vecs-np.eye(len(evals))))
        assert max(residuals)<1e-9 and gram<1e-8
        sectors.append(dict(charge=q,dimension=len(selected),energy=float(evals[0]),energies=evals.tolist(),
            residuals=residuals.tolist(),gram_error=gram,hermiticity_error=0.,runtime_seconds=time.time()-t,
            multiplicity_scope='Returned independently orthogonal eigenpairs; scalar ARPACK does not certify complete multiplicities.'))
        print(family,'Nup',nup,'winding',q,'E',evals[0],'residual',max(residuals),flush=True)
        record=dict(family=family,L=2,width=1,physical_spins=n,nup=nup,total_sz=nup-n/2,filling_fraction=nup/n,
            filling_selection=('Fixed low-number algebra fixture; not a ground filling' if fixture else 'Fixed independently audited finite ground number; not a thermodynamic filling claim'),
            geometry=geometry,loop_label=cycle['label'],edge_exponents=p,physical_order_winding_weights=weights.tolist(),
            full_number_block_dimension=len(basis),sectors=sectors)
        output=dict(status='partial',records=records+[record],runtime_seconds=time.time()-started,
            audit=provenance(seed,'Independent real sparse microscopic winding-block ED',dict(k=8,ncv=80,tolerance=1e-12,residual_gate=1e-9),'Deterministic real random Lanczos start',['Physical N_up','Exact microscopic winding charge']),
            interpretation='Sector-resolved small-cylinder optimizer gate. No original result overwritten; no complete degeneracy, MES, infinite filling or 2D-phase claim.')
        output_path='results/winding_real_sparse_fixture_ed.json' if fixture else 'results/winding_physical_sector_ed.json'
        pathlib.Path(output_path).write_text(json.dumps(output,indent=2)+'\n')
    assert sum(s['dimension'] for s in sectors)==len(basis)
    record['sector_dimensions_cover_full_number_block']=True
    if fixture:
        original=next(r for r in json.load(open('results/sector_penalty_ed_audit.json'))['records'] if r['family']==family)
        assert all(abs(s['energy']-original['sectors'][s['charge']]['energy'])<1e-9 for s in sectors)
    else:
        path='results/zigzag_L2_w1_block_ed.json' if family=='zigzag' else 'results/armchair_L2_w1_Nup8_block_ed.json'
        original=json.load(open(path))
        assert abs(min(s['energy'] for s in sectors)-original['energies'][0])<1e-9
    records.append(record)
output['status']='complete'
output['records']=records
output['runtime_seconds']=time.time()-started
pathlib.Path(output_path).write_text(json.dumps(output,indent=2)+'\n')
