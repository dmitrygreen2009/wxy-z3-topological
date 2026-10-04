"""Lift contractible honeycomb hexagons and export exact microscopic CGS operators."""
import hashlib,json,pathlib
out=pathlib.Path('geometry/cgs_plaquettes');out.mkdir(exist_ok=True)
count=0
for path in sorted(pathlib.Path('geometry').glob('*_L*_w*_*.json')):
    g=json.loads(path.read_text());c=g.get('conventions',{});family=c.get('family')
    if family not in ('zigzag','armchair'):continue
    L=c['length_slices'];w=c['width_cells'];nv=len(g['vertices'])
    lookup={(v['t'],v['j'],v['sublattice']):v['id'] for v in g['vertices']}
    inverse={physical:j+1 for j,physical in enumerate(g['mps_order_physical_spin_ids'])}
    records=[]
    for t in range(1 if family=='zigzag' else 2,L):
        for j in range(w):
            x,y=(t,j) if family=='zigzag' else (t+j,j)
            descriptors=[(x,y,1),(x,y+1,3),(x,y+1,2),(x-1,y+1,1),(x-1,y+1,3),(x,y,2)]
            p=[0]*len(g['shared_gauge_edges']);walk=[]
            for k,(ax,ay,leg) in enumerate(descriptors):
                at,aj=(ax,ay%w) if family=='zigzag' else (ax-ay,ay%w)
                vertex=lookup[(at,aj,'A')];edge=g['endpoint_leg_table'][vertex-1][leg-1]
                direction=1 if k%2==0 else -1;p[edge-1]=(p[edge-1]+direction)%3
                a=[ax,ay,'A'];b=[ax-(leg==2),ay-(leg==3),'B']
                walk.append({'edge_id':edge,'local_leg':leg,'from_unwrapped':a if direction==1 else b,'to_unwrapped':b if direction==1 else a})
            assert all(walk[k]['to_unwrapped']==walk[k+1]['from_unwrapped'] for k in range(5))
            assert walk[0]['from_unwrapped']==walk[-1]['to_unwrapped']
            assert all(sum(p[e-1] for e in row)%3==0 for row in g['endpoint_leg_table'])
            record={'label':f'contractible_hexagon_t{t}_j{j}','edge_exponents':p,'walk':walk,'winding_number':0,'unwrapped_displacement':[0,0],'contractible':True}
            if L==4 and c['ordering']=='star' and t==2:
                leading=g['physical_spins']-36*w;triplets=[];gauges=[]
                for v,row in enumerate(g['endpoint_leg_table'],1):
                    cc=p[row[0]-1];q=(p[row[1]-1]-cc)%3
                    if cc or q:
                        positions=[inverse[3*(v-1)+a]-leading for a in (1,2,3)]
                        assert positions==list(range(positions[0],positions[0]+3))
                        triplets.append([positions[0],cc,q])
                for edge,exponent in enumerate(p,1):
                    if exponent:gauges.append([inverse[3*nv+edge]-leading,exponent])
                record['infinite_star_operator']={'triplets':triplets,'gauges':gauges,'first_site':min(r[0] for r in triplets),'last_site':max([r[0]+2 for r in triplets]+[r[0] for r in gauges]),'coordinate_origin':'First bulk spatial slice t=0; scalar indices translate beyond the two-slice infinite cell'}
            records.append(record);count+=1
    definition={'parent_geometry':str(path),'parent_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'family':family,'L':L,'width':w,'index_base':1,'plaquettes':records}
    (out/path.name).write_text(json.dumps(definition,indent=2)+'\n')
print('Exported and checked',count,'closed contractible hexagon lifts')
