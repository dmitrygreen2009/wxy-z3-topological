"""Exact structural/arithmetic checks; not a numerical phase diagnosis."""
import json
from fractions import Fraction
from pathlib import Path
p=Path(__file__).resolve().parents[1]
g=json.loads((p/'geometry/primitive_honeycomb_cell.json').read_text())
m=[i for v in g['vertices'] for i in v['matter_spins']]
e=g['physical_gauge_edges']
assert sorted(m+[x['spin'] for x in e])==list(range(1,10))
assert len(m)==6 and len(e)==3
assert [x['endpoints'][1]['cell'] for x in e]==[[0,0],[-1,0],[0,-1]]
assert all(x['endpoints'][0]['leg']==x['endpoints'][1]['leg']==i+1 for i,x in enumerate(e))
assert g['B_incident_gauge_owners']==[[0,0],[1,0],[0,1]]
assert Fraction(9,2)%1==Fraction(1,2)
allowed=[Fraction(i,3) for i in range(3)]
assert all((3*q)%1==0 for q in allowed)
assert Fraction(1,2) not in allowed
assert 9*Fraction(4,9)==4 and 9*Fraction(5,9)==5
r={'validation':'passed','primitive_physical_spins':9,'half_filling_per_cell':'9/2','D_Z3_allowed_U1_fractional_charges':[str(q) for q in allowed],'conditional_half_filling_incompatibility':True,'assumptions':['gapped phase','unbroken physical U1','unbroken primitive translations','ordinary nonpermuting background-anyon SET formulation'],'phase_conclusion':'No unconditional conclusion about the microscopic ground state.','numerical_calculations_repeated':False}
(p/'results/primitive_filling_audit.json').write_text(json.dumps(r,indent=2)+'\n')
print('Primitive geometry, filling and order-three charge arithmetic: passed')
