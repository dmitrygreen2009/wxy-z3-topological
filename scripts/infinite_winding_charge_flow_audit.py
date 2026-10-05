"""Exact microscopic charge flow across nonoverlapping infinite winding cells.
This is an algebraic prerequisite, not a validated infinite optimizer.
"""
import collections, json, pathlib, time
from provenance import provenance

started=time.time();records=[];unavailable=[]
for family in ['zigzag','armchair']:
    for w in range(1,5):
        loop_slices=1 if family=='zigzag' else 2;loop_sites=9*w*loop_slices
        weights=[0]*loop_sites
        for j in range(w):
            row=9*j
            matter=[(row+1,1,2),(row+6,1,2)] if family=='zigzag' else [(row+1,1,2),(row+6,1,1),(row+9*w+1,1,1),(row+9*w+6,1,2)]
            gauges=[(row+4,1),(row+5,2)] if family=='zigzag' else [(row+4,1),(row+5,2),(row+9,2),(row+9*w+4,1)]
            for first,c,q in matter:
                for a in range(3):weights[first+a-1]=(c+q*a)%3
            for site,p in gauges:weights[site-1]=p%3
        assert weights[-1]==0
        for slices in [2,3,6]:
            suffix='' if slices==2 else f'_cell{slices}'
            path=f'geometry/infinite_{family}_w{w}{suffix}_star.json'
            if not pathlib.Path(path).exists():
                unavailable.append(path)
                continue
            g=json.load(open(path));assert g['cell_slices']==slices
            checked=0;max_groups=0
            for p in g['pairs']:
                if p['a']!=1:continue
                start=p['m'];i=p['i']-1;other=[a for a in range(3) if a!=i];gauge=p['g']
                variants=[[(start+i,-1),(gauge,1)],[(start+i,1),(start+other[0],-1),(start+other[1],-1),(gauge,1)]]
                for terms in variants:
                    flows=collections.defaultdict(int)
                    for site,dn in terms:
                        group=(site-1)//loop_sites;weight=weights[(site-1)%loop_sites]
                        flows[group]+=weight*dn
                    assert all(value%3==0 for value in flows.values()),(path,p,terms,dict(flows))
                    assert sum(dn for _,dn in terms)==0
                    checked+=1;max_groups=max(max_groups,len(flows))
            records.append(dict(family=family,width=w,source_geometry=path,source_cell_slices=slices,
                nonoverlapping_winding_cell_slices=loop_slices,individual_sites_per_winding_cell=loop_sites,
                onsite_winding_weights_one_loop_cell=weights,checked_distinct_number_flow_terms=checked,
                maximum_winding_cells_touched=max_groups,all_individual_loop_cell_charge_flows_zero=True,
                source_cell_matches_charge_period=(slices%loop_slices==0)))
out=dict(records=records,unavailable_geometry_records=unavailable,runtime_seconds=time.time()-started,
    audit=provenance(None,'Exact termwise microscopic winding-charge algebra',dict(modulus=3,algebra='Validated real eight-state exchange identity'),'Recorded infinite endpoint pairs and explicit microscopic loop actions',['Physical total N_up','Each disjoint microscopic winding generator']),
    conclusion='Every exact exchange product term is neutral separately in each chosen nonoverlapping winding cell. A block constraint at each such boundary is therefore a possible exact optimization subspace. The official infinite update, expansion, canonicalization, and boundary-QN support must still be validated; a globally repeated charge alone is insufficient.',
    caveats=['Armchair odd-slice source cells do not match the two-slice loop period and require a larger cell for this construction.',
        'Fixing a pattern of disjoint winding generators can constrain contractible CGS combinations; it is not automatically one emergent flux or a MES.',
        'No filling, infinite-state convergence, or 2D phase is established.'])
pathlib.Path('results/infinite_winding_charge_flow_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print('Checked',len(records),'geometry/cell combinations: every individual winding-cell charge flow is exactly zero')
