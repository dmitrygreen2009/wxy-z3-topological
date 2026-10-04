"""Separate convergence diagnostics from circumference intercept estimates."""
import glob,json,pathlib
import numpy as np
from provenance import provenance
from fit_entropy import fit,LN3
def bound_state_consistency(r):
    fields=['canonical_error','left_isometry_error','right_isometry_error','center_normalization_error','transfer_normalization_error']
    return all(isinstance(r.get(k),(int,float)) and np.isfinite(r[k]) and r[k]<1e-10 for k in fields)
bond_audit={r["source_result"]:r for r in json.load(open("results/infinite_bond_dimension_audit.json"))["records"]}
points=[]
for path in sorted(glob.glob('results/infinite_*_chi*.json')):
    r=json.load(open(path))
    if 'spatial_entropy' not in r or r.get('validation_fixture',False):continue
    points.append(dict(source_file=path,family=r['family'],width=r['width'],circumference=r.get('physical_circumference',r['width']*(np.sqrt(3) if r['family']=='zigzag' else 3.0)),
        circumference_source='saved measurement' if 'physical_circumference' in r else 'derived from the explicit validated lattice convention',
        cell_slices=r.get('cell_slices',2),cap=r['cap'],historical_reported_chi=r['chi'],actual_chi=bond_audit.get(path,{}).get('maximum_bond_dimension',r.get('bond_dimension_details',{}).get('maximum_bond_dimension')),bond_dimension_audit_status=bond_audit.get(path,{}).get('status','Not yet independently audited'),ordering=r.get('infinite_ordering','matter_first'),
        conserving_ansatz=r.get('u1_conserving_ansatz'),energy_cell=r['energy_cell'],
        number_background=r.get('number_background'),observed_number_density=r.get('mean_sz',0)+.5 if 'mean_sz' in r else None,
        number_density_error_vs_background=r.get('number_density_error_vs_background'),
        disjoint_dimer_variational_upper_bound_cell=-2*r.get('cell_slices',2)*r['width']/np.sqrt(3),
        isolated_A_star_variational_upper_bound_cell=-2.4030921*r.get('cell_slices',2)*r['width'],
        bound_comparison_state_consistency_checked=bound_state_consistency(r),
        trial_bound_comparison_interpretation='Raw above-trial comparison is not an exclusion proof for inconsistent center tensors; fixed-density comparisons also require a compatible trial density.',
        energy_above_known_trial_state=bool(r['energy_cell']>-2.4030921*r.get('cell_slices',2)*r['width']+1e-9),
        entropy=r.get('spatial_entropy_slice_mean',r['spatial_entropy']),entropy_modulation=r.get('slice_cut_entropy_modulation'),
        xi_cells=r['xi_cells'],solver_residual=r.get('solver_residual'),canonical_error=r.get('canonical_error'),
        transfer_residuals=r.get('transfer_residuals'),fixed_point_rank_detected=r.get('fixed_point_rank_detected'),
        transfer_gap_above_residual_scale=r.get('transfer_gap_above_residual_scale'),
        measurement_tag=r.get('measurement_tag',''),
        branch_id=r.get('measurement_tag','').removesuffix('_fixedpoint'),
        sector_consistency_certified=False,physical_bulk_gap_certified=False,
        fixed_point_detection_scope='Two independent starts detect some repeated fixed points; this is not a certified complete peripheral spectrum.',
        geometry_limitations='Narrow parallel-edge quotient has extra local S3 symmetry; inspect cycle-charge state selection' if r['family']=='zigzag' and r['width']==1 else 'Circumference convergence remains required'))
# Preserve every raw measurement. Use fixed-point remeasurement in preference
# to the older center tensors of the same optimized ansatz when available.
selected=[]
keys=set((p['family'],p['width'],p['cap'],p['ordering'],p['conserving_ansatz'],p['branch_id'],p['cell_slices']) for p in points)
for key in sorted(keys,key=str):
    group=[p for p in points if (p['family'],p['width'],p['cap'],p['ordering'],p['conserving_ansatz'],p['branch_id'],p['cell_slices'])==key]
    fixed=[p for p in group if p['measurement_tag'].endswith('_fixedpoint')]
    selected.append(min(fixed or group,key=lambda p:p['energy_cell']))
convergence=[]
for key in sorted(set((p['family'],p['width'],p['ordering'],p['conserving_ansatz'],p['branch_id'],p['cell_slices']) for p in selected),key=str):
    group=sorted([p for p in selected if (p['family'],p['width'],p['ordering'],p['conserving_ansatz'],p['branch_id'],p['cell_slices'])==key],key=lambda p:p['cap'])
    for a,b in zip(group,group[1:]):
        valid=all(isinstance(p['xi_cells'],(int,float)) and np.isfinite(p['xi_cells']) and p['xi_cells']>0 for p in [a,b])
        log_ratio=float(np.log(b['xi_cells']/a['xi_cells'])) if valid else None
        convergence.append(dict(family=key[0],width=key[1],ordering=key[2],conserving_ansatz=key[3],
            branch_id=key[4],cell_slices=key[5],sector_consistency_certified=False,lower=a,higher=b,entropy_change=b['entropy']-a['entropy'],energy_change=b['energy_cell']-a['energy_cell'],
            xi_ratio=b['xi_cells']/a['xi_cells'] if valid else None,
            finite_entanglement_effective_c=6*(b['entropy']-a['entropy'])/log_ratio if valid and abs(log_ratio)>1e-12 else None,
            interpretation='Two-point finite-entanglement slope only; unconverged solvers or a narrow circumference cannot establish a critical bulk phase.'))
fits=[]
for key in sorted(set((p['family'],p['cap'],p['ordering'],p['conserving_ansatz'],p['branch_id'],p['cell_slices']) for p in selected),key=str):
    group=sorted([p for p in selected if (p['family'],p['cap'],p['ordering'],p['conserving_ansatz'],p['branch_id'],p['cell_slices'])==key and not p['energy_above_known_trial_state']],key=lambda p:p['width'])
    if len(group)<3:continue
    result=fit(group)
    result.update(family=key[0],cap=key[1],ordering=key[2],conserving_ansatz=key[3],branch_id=key[4],cell_slices=key[5],raw_points=group,eligible_for_topological_inference=False,
        smallest_circumference_removed=fit(group[1:]),status='Exploratory until solver, bond dimension, state sector and circumference converge.')
    result['smallest_removed_gamma_change']=result['smallest_circumference_removed']['gamma']-result['gamma']
    fits.append(result)
output=dict(ln3=LN3,raw_measurements=points,selected_measurements=selected,bond_convergence=convergence,circumference_fits=fits,accepted_topological_fits=[],
    scientific_strategy_reaudit='results/scientific_strategy_reaudit.md',
    audit=provenance(None,'Saved infinite-MPS entropy/transfer analysis',{},'Recorded full-state infinite measurements',[]),
    status='No converged topological intercept is inferred automatically. Empty fits mean insufficient distinct circumferences.')
pathlib.Path('results/infinite_analysis.json').write_text(json.dumps(output,indent=2))
print('Infinite measurements:',len(points),'convergence comparisons:',len(convergence),'circumference fits:',len(fits))
