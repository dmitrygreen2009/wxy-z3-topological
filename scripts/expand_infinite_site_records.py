"""Export explicit physical-site and endpoint tables from checked pair records."""
import json,pathlib,math

def decode(index,w,ordering):
 t,k=divmod(index-1,9*w);k+=1
 if ordering=='star':
  j,r=divmod(k-1,9);r+=1
  if r<=3:return dict(kind='matter',t=t,j=j,sublattice='A',species=r)
  if 6<=r<=8:return dict(kind='matter',t=t,j=j,sublattice='B',species=r-5)
  return dict(kind='gauge',owner_t=t+(r==9),owner_j=j,leg={4:1,5:3,9:2}[r])
 if k<=6*w:
  j,r=divmod(k-1,6);return dict(kind='matter',t=t,j=j,sublattice='A' if r<3 else 'B',species=r%3+1)
 j,r=divmod(k-6*w-1,3);return dict(kind='gauge',owner_t=t+(r==2),owner_j=j,leg=[1,3,2][r])

for path in sorted(pathlib.Path('geometry').glob('infinite_*_w*.json')):
 r=json.loads(path.read_text());w=r['width'];n=r['cell_spins'];family=r['family'];ordering=r['ordering']
 sites=[dict(index=i,**decode(i,w,ordering)) for i in range(1,n+1)]
 vertices={};edges=[]
 for s in sites:
  if s['kind']=='matter':
   key=(s['t'],s['j'],s['sublattice']);v=vertices.setdefault(key,dict(t=key[0],j=key[1],sublattice=key[2],matter_sites=[],incident_gauge_sites=[]));v['matter_sites'].append(s['index'])
  else:
   t,j,i=s['owner_t'],s['owner_j'],s['leg'];bt=t-(i==2)+(i==3 and family=='armchair');bj=(j-(i==3))%w
   edges.append(dict(physical_site=s['index'],owner_A=[t,j],local_leg=i,endpoints=[dict(sublattice='A',t=t,j=j,leg=i),dict(sublattice='B',t=bt,j=bj,leg=i)]))
 for p in r['pairs']:
  if p['a']!=1:continue
  shift=((p['m']-1)//n)*n;m=decode(p['m']-shift,w,ordering);g=p['g']-shift
  vertices[(m['t'],m['j'],m['sublattice'])]['incident_gauge_sites'].append(dict(leg=p['i'],site_in_cell=(g-1)%n+1,cell_translation=(g-1)//n))
 assert len(vertices)==2*w*r['cell_slices'] and len(edges)==3*w*r['cell_slices']
 assert all(len(v['matter_sites'])==3 and sorted(x['leg'] for x in v['incident_gauge_sites'])==[1,2,3] for v in vertices.values())
 r.update(cell_translation_primitive_xy=[r['cell_slices'],0],cell_translation_physical_xy=[math.sqrt(3)/2*r['cell_slices'],1.5*r['cell_slices']],cell_translation_axial_projection=r['cell_slices']*(1.5 if family=='zigzag' else math.sqrt(3)/2),cell_translation_circumferential_projection=r['cell_slices']*(math.sqrt(3)/2 if family=='zigzag' else 1.5),translation_interpretation='Primitive-a1 cell translation includes a circumferential shift; xi_physical_axial uses its perpendicular projection, not the Euclidean vector length.',physical_sites=sites,vertices=list(vertices.values()),physical_shared_gauge_edges=edges,periodic_identifications=dict(circumference_primitive_translation=[0,w] if family=='zigzag' else [w,w],axial_slice_coordinate='x' if family=='zigzag' else 'x-y',cell_translation_slices=r['cell_slices']))
 path.write_text(json.dumps(r,indent=2)+'\n')
print('Expanded physical-site, vertex and shared-edge records for all saved infinite cells')
