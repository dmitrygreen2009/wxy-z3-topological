import json,pathlib
for path in pathlib.Path('geometry').glob('*.json'):
    r=json.load(open(path))
    if 'shared_gauge_edges' not in r:continue
    nv=len(r['vertices']);matter=r['matter_spins'];edges=r['shared_gauge_edges']
    assert len(matter)==3*nv
    ids=[m['physical_spin_id'] for m in matter]+[e['physical_spin_id'] for e in edges]
    assert sorted(ids)==list(range(1,r['physical_spins']+1))
    for e in edges:
        assert len(e['endpoints']) in (1,2)
        for p in e['endpoints']:
            assert r['endpoint_leg_table'][p['vertex_id']-1][p['local_leg_i']-1]==e['edge_id']
    assert sum(len(e['endpoints']) for e in edges)==3*nv
for name,n in [('single_star',6),('shared_edge_pair',11),('torus_four_star',18)]:
    assert json.load(open(f'geometry/{name}.json'))['physical_spins']==n
r=json.load(open('geometry/torus_four_star.json'))
pairs=[tuple(sorted(p['vertex_id'] for p in e['endpoints'])) for e in r['shared_gauge_edges']]
assert len(set(pairs))==4 and len(pairs)==6
print('Explicit geometry records passed, including distinct parallel torus edges.')
