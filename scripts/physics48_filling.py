"""Exact rational fusion/translation checks; no many-body calculation."""
from fractions import Fraction
import itertools,json,pathlib
import numpy as np
anyons=list(itertools.product(range(3),repeat=2))
autos=[]
for entries in itertools.product(range(3),repeat=4):
 M=np.array(entries).reshape(2,2)
 if round(np.linalg.det(M))%3==0:continue
 if all(int(np.prod(M@a%3))%3==(a[0]*a[1])%3 for a in anyons):autos.append(M)
assert len(autos)==4
records=[]
for i,tx in enumerate(autos):
 for j,ty in enumerate(autos):
  assert np.array_equal(tx@ty%3,ty@tx%3)
  charges=[]
  for qa,qm in anyons:
   Q=lambda a:Fraction((qa*a[0]+qm*a[1])%3,3)
   if all(Q(tx@a%3)==Q(a) and Q(ty@a%3)==Q(a) for a in anyons):
    assert all(3*Q(a)%1==0 for a in anyons)
    assert all(Q(a)!=Fraction(1,2) for a in anyons)
    charges.append([qa,qm])
  records.append(dict(tx=tx.tolist(),ty=ty.tolist(),translation_invariant_charge_characters=charges))
primitive=json.load(open('results/primitive_filling_audit.json'));bound=json.load(open('results/number_bound_audit.json'))
out=dict(primitive_spins=9,half_filling_per_primitive_cell='9/2',local_physical_charge='S^z+1/2 is integer-valued 0/1',
 allowed_anyon_charges=['0','1/3','2/3'],braided_translation_permutation_cases=records,
 conditional_no_go_including_anyonic_translation_permutations=True,
 assumptions=['gapped ordinary bosonic D(Z3) SET', 'unbroken onsite physical U1', 'unbroken primitive translations', 'half filling', 'no additional symmetry-breaking ground degeneracy', 'ordinary commuting translations, no external magnetic translation algebra'],
 proof='Connected U1 does not permute anyons. Its 2pi flux anyon v lies in A=Z3xZ3 and obeys v^3=1. Translation preserves charge, hence v. Three flux insertions are scalar on the purely topological ground space, but Oshikawa translation algebra gives exp(2pi i*3*nu*Ly). At nu=9/2 and odd Ly this is -1, a contradiction. Equivalently the required fractional-charge half spinon is absent: all anyons have 3Qa integer.',
 flux_argument_status='Conditional no-go under standard gapped SET flux-threading assumptions; does not prove the actual model has a gap, half filling, or unbroken translations',
 previous_audit_scope='Prior nonpermuting audit retained; this extraction checks continuous-U1 permutations explicitly',
 number_bound_window=bound['bulk_density_window_from_exact_lower_bound_and_trial'],
 allowed_density_grid=primitive['ordinary_D_Z3_symmetric_gapped_density_grid'],
 alternatives=['different filling', 'primitive translation breaking (even cell enlargement removes half anomaly)', 'U1 breaking', 'gaplessness', 'topological order with appropriate even-denominator charge structure or additional order'],
 source_links=['https://arxiv.org/abs/1511.02263','https://arxiv.org/abs/1410.2894','https://arxiv.org/abs/1808.00324'],many_body_calculations_repeated=False)
pathlib.Path('results/physics48_filling.json').write_text(json.dumps(out,indent=2)+'\n')
print('All 16 braided translation-permutation pairs checked; none supplies half U1 charge.')
