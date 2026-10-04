"""Compare saved full-state density profiles in physical axial coordinates."""
import argparse,collections,hashlib,json,pathlib,subprocess
p=argparse.ArgumentParser();p.add_argument('neutral');p.add_argument('charged');p.add_argument('output');a=p.parse_args()
n=json.loads(pathlib.Path(a.neutral).read_text());c=json.loads(pathlib.Path(a.charged).read_text())
assert (n['family'],n['length'],n['width'],n['ordering'])==(c['family'],c['length'],c['width'],c['ordering'])
f=n['family'];L=n['length'];w=n['width']
g=json.loads(pathlib.Path(f'geometry/{f}_L{L}_w{w}_{n["ordering"]}.json').read_text())
vertices={v['id']:v for v in g['vertices']};u={}
for m in g['matter_spins']:
 v=vertices[m['vertex_id']];u[m['physical_spin_id']]=v['t']+(1/3 if f=='zigzag' and v['sublattice']=='B' else 0)
for edge in g['shared_gauge_edges']:
 x,y=edge['owner_A_unwrapped_bravais'];at=x if f=='zigzag' else x-y
 leg=edge['endpoints'][0]['local_leg_i'];u[edge['physical_spin_id']]=at+({1:1/6,2:-1/3,3:1/6} if f=='zigzag' else {1:0,2:-.5,3:.5})[leg]
nd=dict(zip(n['mps_order_physical_spin_ids'],n['site_sz_profile_mps_order']));cd=dict(zip(c['mps_order_physical_spin_ids'],c['site_sz_profile_mps_order']))
assert set(nd)==set(cd)==set(u)
rows=[{'physical_spin_id':j,'physical_axial_slice_coordinate':u[j],'neutral_sz':nd[j],'charged_sz':cd[j],'difference_charged_minus_neutral':cd[j]-nd[j]} for j in sorted(u)]
difference=sum(r['difference_charged_minus_neutral'] for r in rows)
assert abs(difference-(c['actual_nup']-n['actual_nup']))<1e-8
lo,hi=L/4,3*L/4
central=sum(r['difference_charged_minus_neutral'] for r in rows if lo<=u[r['physical_spin_id']]<hi)
result={'family':f,'length':L,'width':w,'neutral_source':a.neutral,'charged_source':a.charged,
 'source_sha256':{a.neutral:hashlib.sha256(pathlib.Path(a.neutral).read_bytes()).hexdigest(),a.charged:hashlib.sha256(pathlib.Path(a.charged).read_bytes()).hexdigest()},
 'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'density_difference_sum':difference,'central_half_open_interval':[lo,hi],
 'central_region_number_difference':central,'complement_number_difference':difference-central,'physical_spin_profiles':rows,
 'interpretation':'Finite-state density differences and region sums only. Number and loop sectors can both differ. A short cylinder cannot certify whether a thermodynamic filling shift is a boundary or bulk effect; length and circumference convergence are still required.'}
pathlib.Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print('total delta N:',difference,'central:',central,'complement:',difference-central)
