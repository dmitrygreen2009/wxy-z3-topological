"""Static independent audit of embedded cuts, infinite incidence and cycle winding.
No ED, DMRG or VUMPS optimization is repeated.
"""
import collections,json,pathlib,hashlib
from provenance import provenance
reports=[];cycle_dir=pathlib.Path('geometry/cgs_cycles');cycle_dir.mkdir(exist_ok=True)
for path in sorted(pathlib.Path('geometry').glob('*_L*_w*_*.json')):
    graph=json.loads(path.read_text());c=graph.get('conventions',{})
    family=c.get('family')
    if family not in ('zigzag','armchair'):continue
    L=c['length_slices'];w=c['width_cells'];nv=len(graph['vertices']);edges=graph['shared_gauge_edges']
    expected_matter=[0,1/3] if family=='zigzag' else [0,0]
    expected_gauge=[1/6,-1/3,1/6] if family=='zigzag' else [0,-.5,.5]
    assert c['physical_axial_matter_offsets_A_B']==expected_matter
    assert c['physical_axial_gauge_offsets_from_owner_A_by_leg']==expected_gauge
    assert nv==2*L*w and len(edges)==3*L*w+(w if family=='zigzag' else 2*w)
    vertices={v['id']:v for v in graph['vertices']};expected=set()
    for matter in graph['matter_spins']:
        v=vertices[matter['vertex_id']];u6=6*v['t']+(2 if family=='zigzag' and v['sublattice']=='B' else 0)
        if u6<3*L:expected.add(matter['physical_spin_id'])
    for edge in edges:
        ax,ay=edge['owner_A_unwrapped_bravais'];t=ax if family=='zigzag' else ax-ay
        legs={p['local_leg_i'] for p in edge['endpoints']};assert len(legs)==1
        leg=next(iter(legs));offset={1:1,2:-2,3:1} if family=='zigzag' else {1:0,2:-3,3:3}
        if 6*t+offset[leg]<3*L:expected.add(edge['physical_spin_id'])
        if 'periodic_translation_to_B_representative' in edge:
            dx,dy=edge['periodic_translation_to_B_representative']
            assert (dx==0 and dy%w==0) if family=='zigzag' else (dx==dy and dx%w==0)
    actual=set(graph['mps_order_physical_spin_ids'][:graph['spatial_cut_mps_bond']]);assert expected==actual
    lookup={(v['t'],v['j'],v['sublattice']):v['id'] for v in vertices.values()};cycles=[]
    # Lift the walk to the unwrapped honeycomb. Its displacement, rather than
    # the number of quotient edges, determines contractibility.
    for t in range(L if family=='zigzag' else L-1):
        p=[0]*len(edges);walk=[]
        for j in range(w):
            descriptors=[(t,j,1,1),(t,j+1,3,-1)] if family=='zigzag' else [(t,j,1,1),(t+1,j,2,-1),(t+1,j,1,1),(t,j+1,3,-1)]
            for at,aj,leg,direction in descriptors:
                vertex=lookup[(at,aj%w,'A')];e=graph['endpoint_leg_table'][vertex-1][leg-1]
                p[e-1]=(p[e-1]+direction)%3
                x,y=(at,aj) if family=='zigzag' else (at+aj,aj)
                a=[x,y,'A'];b=[x-(leg==2),y-(leg==3),'B']
                walk.append(dict(edge_id=e,local_leg=leg,from_unwrapped=a if direction==1 else b,to_unwrapped=b if direction==1 else a))
        assert all(walk[k]['to_unwrapped']==walk[k+1]['from_unwrapped'] for k in range(len(walk)-1))
        start=walk[0]['from_unwrapped'];end=walk[-1]['to_unwrapped'];displacement=[end[k]-start[k] for k in [0,1]]
        assert displacement==([0,w] if family=='zigzag' else [w,w])
        assert all(sum(p[e-1] for e in row)%3==0 for row in graph['endpoint_leg_table'])
        cycles.append(dict(label=f'circumference_cycle_t{t}',edge_exponents=p,winding_number=1,
            unwrapped_displacement=displacement,walk=walk,contractible=False,
            interpretation='Exact microscopic CGS cycle; its eigencharge is not by itself a proof of topological order or of a minimally entangled ground state.'))
    definition=dict(parent_geometry=str(path),parent_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),parent_semantic_sha256=hashlib.sha256(json.dumps(graph,sort_keys=True,separators=(',',':')).encode()).hexdigest(),family=family,L=L,width=w,
        index_base=1,cycles=cycles,narrow_zigzag_note='At width one this noncontractible circumference cycle is a two-edge parallel cycle; it has extra local S3 symmetry.' if family=='zigzag' and w==1 else None)
    (cycle_dir/path.name).write_text(json.dumps(definition,indent=2)+'\n')
    reports.append(dict(file=str(path),family=family,L=L,width=w,physical_spins=graph['physical_spins'],physical_cut_passed=True,cycle_count=len(cycles)))

def decode(index,w,ordering):
    t,k=divmod(index-1,9*w);k+=1
    if ordering=='star':
        j,r=divmod(k-1,9);r+=1
        if r<=3:return ('matter',t,j,'A',r)
        if 6<=r<=8:return ('matter',t,j,'B',r-5)
        return ('gauge',t+(r==9),j,{4:1,5:3,9:2}[r])
    if k<=6*w:
        j,r=divmod(k-1,6);return ('matter',t,j,'A' if r<3 else 'B',r%3+1)
    j,r=divmod(k-6*w-1,3);return ('gauge',t+(r==2),j,[1,3,2][r])
infinite=[]
for path in sorted(pathlib.Path('geometry').glob('infinite_*_w*.json')):
    table=json.loads(path.read_text());family=table['family'];w=table['width'];ordering=table['ordering'];n=18*w
    assert table['cell_spins']==n and len(table['pairs'])==36*w
    signatures=collections.Counter()
    for pair in table['pairs']:
        m=decode(pair['m'],w,ordering);g=decode(pair['g'],w,ordering)
        assert m[0]=='matter' and g[0]=='gauge' and m[-1]==pair['a'] and g[-1]==pair['i']
        _,t,j,s,a=m;x,y=(t,j) if family=='zigzag' else (t+j,j)
        if s=='B':x+=pair['i']==2;y+=pair['i']==3
        at,aj=(x,y%w) if family=='zigzag' else (x-y,y%w)
        assert g[1:]==(at,aj,pair['i'])
        signatures[(t%2,j,s,a,pair['i'])]+=1
    assert len(signatures)==36*w and set(signatures.values())=={1}
    for boundary in [9*w,18*w]:
        for index in range(-n+1,2*n+1):
            spin=decode(index,w,ordering)
            if spin[0]=='matter':u6=6*spin[1]+(2 if family=='zigzag' and spin[3]=='B' else 0)
            else:u6=6*spin[1]+({1:1,2:-2,3:1} if family=='zigzag' else {1:0,2:-3,3:3})[spin[3]]
            assert (index<=boundary)==(u6<6*(boundary//(9*w)))
    infinite.append(dict(file=str(path),family=family,width=w,ordering=ordering,incidence_passed=True,spatial_slice_cuts_passed=True))
output=dict(finite_geometry_checks=reports,infinite_geometry_checks=infinite,
    interpretation='Static coordinate, edge-identification, spatial-cut and winding audit. No validated numerical calculation was rerun.',
    audit=provenance(None,'Independent embedded-geometry and homology audit',{},'Committed machine-readable physical lattice definitions',[]))
pathlib.Path('results/strategy_geometry_audit.json').write_text(json.dumps(output,indent=2)+'\n')
print('Physical-cut checks:',len(reports),'infinite incidence/cut checks:',len(infinite),'circumference cycles exported:',sum(r['cycle_count'] for r in reports))
