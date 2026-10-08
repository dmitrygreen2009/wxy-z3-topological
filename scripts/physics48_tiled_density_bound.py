"""Combine a measured valid 2D trial with already validated operator inequalities."""
import json,math,pathlib
from fractions import Fraction
r=json.load(open('results/physics48_tiled_trial.json'));b=json.load(open('results/physics48_density_bound.json'))
assert r['patch_physical_spin_count']==9*r['primitive_cells_per_patch']
# Explicit upward numerical safety allowance, never a changed solver tolerance.
margin=1e-8;upper=math.ceil((r['bulk_trial_energy_per_primitive_cell']+margin)*1e8)/1e8
lo=-upper/(9*math.sqrt(2));clo=(upper-2*b['cluster_epsilon_conservative'])/(9*b['cluster_multiplier'])
assert lo>b['analytic_density_window'][0] and clo>b['cluster_density_window'][0]
out=dict(measured_2D_trial_energy_per_primitive_cell=r['bulk_trial_energy_per_primitive_cell'],upward_rounded_trial_upper=upper,trial_safety_allowance=margin,number_operator_density_window=[lo,1-lo],same_local_cluster_density_window=[clo,1-clo],local_cluster_multiplier_reused=b['cluster_multiplier'],local_cluster_epsilon_reused=b['cluster_epsilon_conservative'],ordinary_symmetric_D_Z3_density_grid_remaining=[str(Fraction(k,27)) for k in range(28) if clo<=k/27<=1-clo],state_phase_claim=False,benchmark_and_measurement_sources=['results/physics48_tiled_trial.json','results/physics48_density_bound.json','results/physics48_tiled_trial_validation.json'],qualification='Exact physical tiling construction and operator inequalities; numerical endpoint evaluation uses stated margins and independent Bell-state fixtures, not interval-arithmetic certification. No shifted-star eigensolve repeated.')
pair=json.load(open('results/benchmark_number_sectors.json'))['shared_edge_pair']
assert pair['complete_number_scan'] and max(x['residual'] for x in pair['points'])<1e-8
pair_energy=min(x['energy'] for x in pair['points']);lower=math.floor((pair_energy-margin)*1e8)/1e8
out['primitive_cell_pair_cover_lower_bound']=dict(benchmark_minimum=pair_energy,downward_rounded_bound=lower,ground_energy_per_primitive_cell_window=[lower,upper],proof='Pair each primitive A and B across leg0 in the full 2D honeycomb. H is the sum of these overlapping eleven-spin pair operators, each bounded below by the preserved full-number pair minimum. External gauges are one physical spin shared between pair blocks. This graph decomposition applies to the full 2D lattice, not a tiny quotient with extra parallel edges.')
pathlib.Path('results/physics48_tiled_density_bound.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
