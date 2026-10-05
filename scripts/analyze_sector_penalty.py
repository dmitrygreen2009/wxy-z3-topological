"""Compare converged fixed-number microscopic winding-sector MPS candidates.
No entropy point is admitted to a TEE fit by this validation alone.
"""
import csv,glob,json,pathlib
ed={}
for path in ['results/sector_penalty_ed_audit.json','results/winding_physical_sector_ed.json']:
    if pathlib.Path(path).exists():
        ed.update({(r['family'],r['L'],r['width'],r['nup'],s['charge']):s['energy']
            for r in json.load(open(path))['records'] for s in r['sectors']})
rows=[]
for path in sorted(glob.glob('results/*_winding*_penalty_seed*.json')+glob.glob('results/*_winding*_qn_seed*.json')):
    d=json.load(open(path));r=d['records'][-1];q=d['audit']['solver_settings']['winding_charge']
    row=dict(source=path,algorithm='direct winding QN' if d.get('basis')=='exact_matter_charge_basis' else 'exact projector penalty',family=d['family'],L=d['length'],width=d['width'],physical_spins=d['spins'],
        nup=d['nup'],total_sz=d['total_sz'],filling_fraction=d['filling_fraction'],filling_selection=d['filling_selection'],
        winding_charge=q,loop_label=d['loop_label'],cap=r['cap'],actual_bond_dimension=r['maxlinkdim'],
        energy=r['energy'],entropy=r['entropy'],energy_variance=r['energy_variance'],
        loop_mean_real=r['loop_expectation_real'],loop_mean_imag=r['loop_expectation_imag'],
        loop_variance=r['loop_variance'],loop_purity_error=r['sector_purity_error'],
        sector_convergence_passed=r['converged_within_fixed_number_and_loop_sector'],
        correlation_length=None,runtime_seconds=d['runtime_seconds'],git_commit=d['git_commit'],
        eligible_for_topological_entropy_fit=False)
    target=ed.get((row['family'],row['L'],row['width'],row['nup'],q))
    row['independent_ed_energy']=target;row['energy_error_vs_ed']=None if target is None else row['energy']-target
    if row['sector_convergence_passed']:
        assert row['loop_purity_error']<1e-10 and row['energy_variance']<1e-8
        assert target is None or abs(row['energy_error_vs_ed'])<1e-9,'MPS disagrees with independent sector ED'
    rows.append(row)
comparisons=[]
for key in sorted(set((r['family'],r['L'],r['width'],r['nup'],r['loop_label'],r['algorithm']) for r in rows)):
    states={r['winding_charge']:r for r in rows if (r['family'],r['L'],r['width'],r['nup'],r['loop_label'],r['algorithm'])==key and r['sector_convergence_passed']}
    if 0 not in states:continue
    for q in [1,2]:
        if q in states:
            comparisons.append(dict(family=key[0],L=key[1],width=key[2],nup=key[3],loop_label=key[4],algorithm=key[5],nontrivial_charge=q,
                energy_difference_nontrivial_minus_trivial=states[q]['energy']-states[0]['energy'],
                entropy_difference_nontrivial_minus_trivial=states[q]['entropy']-states[0]['entropy'],
                sources=[states[0]['source'],states[q]['source']],
                interpretation='Finite-cylinder fixed-number winding-sector comparison; not a deconfinement, MES or 2D topological certificate.'))
output=dict(states=rows,comparisons=comparisons,accepted_topological_entropy_points=[],
    interpretation='Sector constraint validation precedes physical filling, length/circumference and MES convergence. No gamma is fitted here.')
pathlib.Path('results/sector_penalty_comparison.json').write_text(json.dumps(output,indent=2)+'\n')
if rows:
    with open('results/sector_penalty_states.csv','w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
print('Sector candidates:',len(rows),'converged trivial/nontrivial comparisons:',len(comparisons))
