"""Annotate legacy cylinder manifests without changing any physical graph or cut."""
import json,pathlib
items=[]
for p in pathlib.Path('geometry').glob('*_L*_w*_*.json'):
    g=json.loads(p.read_text());c=g.get('conventions',{});f=c.get('family')
    if f not in ('zigzag','armchair'):continue
    assert c['B_minus_A']==[0,1] and c['B_endpoint_uses_same_leg']
    assert c['A_leg_displacements_to_B']==[[0,0],[-1,0],[0,-1]]
    assert c['periodic_translation_xy']==([0,c['width_cells']] if f=='zigzag' else [c['width_cells']]*2)
    items.append((p,g,c,f))
for p,g,c,f in items:
    c['spatial_partition_rule']='Half-open physical axial coordinate u < L/2, in spatial-slice units'
    c['physical_axial_matter_offsets_A_B']=[0,1/3] if f=='zigzag' else [0,0]
    c['physical_axial_gauge_offsets_from_owner_A_by_leg']=[1/6,-1/3,1/6] if f=='zigzag' else [0,-.5,.5]
    c['ordering_partition_proxy']='Matter key t+0.01a; gauge key (owner_A_t+target_B_t)/2+0.1. For integer L these keys give the same half-open partition as the physical coordinates; they are not physical length coordinates.'
    text=json.dumps(g,separators=(',',':'))
    if text!=p.read_text():p.write_text(text)
print('Physical-coordinate metadata checked for',len(items),'cylinder manifests; run the independent geometry audit next.')
