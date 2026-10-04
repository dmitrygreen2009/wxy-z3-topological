"""Measured full-state Schmidt spectra; numerical clustering is descriptive."""
import glob,json,pathlib,re
import numpy as np
from provenance import provenance

def summarize(probabilities):
    p=np.sort(np.asarray(probabilities,dtype=float))[::-1]
    assert np.all(p>=-1e-14) and abs(sum(p)-1)<1e-8
    p=p[p>0];p/=sum(p)
    levels=-np.log(p)
    clusters={}
    for tolerance in [1e-6,1e-4]:
        groups=[];start=0
        for stop in range(1,len(levels)+1):
            if stop==len(levels) or levels[stop]-levels[start]>tolerance:
                groups.append(dict(level=float(np.mean(levels[start:stop])),multiplicity=stop-start,weight=float(sum(p[start:stop]))));start=stop
        clusters[str(tolerance)]=groups[:30]
    return dict(entropy=float(-sum(p*np.log(p))),renyi2=float(-np.log(sum(p*p))),
        schmidt_rank=len(p),largest_probability=float(p[0]),
        leading_probabilities=p[:30].tolist(),leading_entanglement_levels=levels[:30].tolist(),
        numerical_clusters=clusters,weight_outside_leading30=float(sum(p[30:])))

states=[]
for path in sorted(glob.glob('results/*_L*_w*_chi*.json')):
    r=json.load(open(path))
    if 'family' not in r or 'length' not in r:continue
    for k,q in enumerate(r.get('records',[])):
        if 'schmidt_probabilities' not in q:continue
        value=summarize(q['schmidt_probabilities'])
        assert abs(value['entropy']-q['entropy'])<1e-8
        states.append(dict(source_file=path,record_index=k,family=r['family'],length=r['length'],width=r['width'],
            initialization_branch=re.sub(r'_chi\d+','',pathlib.Path(path).stem).removesuffix('_refined'),seed=r.get('seed'),
            circumference=r['physical_circumference'],nup=r.get('nup',(r['spins']+1)//2),number_sector_scope='Global minimum requires independent number-sector search',ordering=r.get('ordering','axial'),cap=q['cap'],
            bond_dimension=q['maxlinkdim'],energy=q['energy'],spectrum=value))
output=dict(states=states,audit=provenance(None,'Analysis of saved full-state Schmidt probabilities',{},'Saved microscopic MPS spatial-cut spectra',[]),
    interpretation='Entanglement levels and approximate degeneracies depend on state, cut, sector and convergence. No universal degeneracy pattern or topological order is assumed.')
pathlib.Path('results/entanglement_spectra.json').write_text(json.dumps(output,indent=2))
print('Validated and summarized',len(states),'full-state Schmidt spectra.')
