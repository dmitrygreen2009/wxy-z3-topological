import copy,json,pathlib,subprocess,sys,tempfile
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
    c=r.get('conventions',{})
    if c.get('family') in ('zigzag','armchair'):
        # Require the exporter schema; absent metadata must not silently default.
        zigzag=c['family']=='zigzag'
        assert c['physical_axial_matter_offsets_A_B']==([0,1/3] if zigzag else [0,0]),path
        assert c['physical_axial_gauge_offsets_from_owner_A_by_leg']==([1/6,-1/3,1/6] if zigzag else [0,-.5,.5]),path
for name,n in [('single_star',6),('shared_edge_pair',11),('torus_four_star',18)]:
    assert json.load(open(f'geometry/{name}.json'))['physical_spins']==n
r=json.load(open('geometry/torus_four_star.json'))
pairs=[tuple(sorted(p['vertex_id'] for p in e['endpoints'])) for e in r['shared_gauge_edges']]
assert len(set(pairs))==4 and len(pairs)==6
audit=pathlib.Path('scripts/strategy_geometry_audit.py').resolve()
fixture=json.loads(pathlib.Path('geometry/armchair_L6_w3_star.json').read_text())
for fault in ('missing_offset','incorrect_offset','incorrect_cut'):
    with tempfile.TemporaryDirectory() as directory:
        root=pathlib.Path(directory);(root/'geometry').mkdir()
        graph=copy.deepcopy(fixture)
        if fault=='missing_offset':
            del graph['conventions']['physical_axial_matter_offsets_A_B']
        elif fault=='incorrect_offset':
            graph['conventions']['physical_axial_matter_offsets_A_B']=[0,1/3]
        else:
            graph['spatial_cut_mps_bond']+=1
        (root/'geometry/armchair_L6_w3_star.json').write_text(json.dumps(graph))
        result=subprocess.run([sys.executable,str(audit)],cwd=root,capture_output=True,text=True)
        assert result.returncode!=0 and 'AssertionError' in result.stderr,(fault,result.stdout,result.stderr)
print('Explicit geometry records passed, including distinct parallel torus edges.')
