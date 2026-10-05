"""Describe optimization gates separately from the evidence needed for a 2D phase."""
import json,math,pathlib,subprocess
from finite_results import finite_groups
from finite_energy_quality import energy_quality_gates
thresholds=dict(last_sweep_energy_change=1e-8,entropy_refinement_change=1e-4,
                truncation_error=1e-8,bond_entropy_change=1e-3,
                bond_energy_per_vertex_change=1e-6,canonical_solver_residual=1e-7,
                measured_canonical_consistency=1e-10,energy_variance=1e-8,known_sector_energy=1e-9)
def small(value,tolerance):
    return None if value is None else bool(math.isfinite(value) and abs(value)<tolerance)
finite=[]
for path,r,indexed in finite_groups():
    stages=[q for _,q in indexed];index,q=indexed[-1];energies=q.get('sweep_energies',[])
    prior=[s for s in stages if s['cap']<q['cap']]
    truncations=q.get('sweep_max_truncation_errors',q.get('truncation_errors',[])) or []
    errors=dict(last_sweep_energy_change=energies[-1]-energies[-2] if len(energies)>1 else None,
                entropy_refinement_change=q.get('entropy_change_refinement'),
                truncation_error=truncations[-1] if truncations else None,
                bond_entropy_change=q['entropy']-prior[-1]['entropy'] if prior else None,
                bond_energy_per_vertex_change=(q['energy']-prior[-1]['energy'])/(2*r['length']*r['width']) if prior else None)
    gates={k:small(v,thresholds[k]) for k,v in errors.items()}
    legacy_passed=all(v is True for v in gates.values())
    quality=energy_quality_gates(q)
    gates['energy_variance']=quality['energy_variance']
    if q.get('known_fixed_sector_ed_energy') is not None:
        gates['known_sector_energy']=quality['known_sector_energy']
    errors['energy_variance']=q.get('energy_variance')
    errors['known_sector_energy']=None if q.get('known_fixed_sector_ed_energy') is None else q['energy']-q['known_fixed_sector_ed_energy']
    finite.append(dict(source_file=path,record_index=index,family=r['family'],length=r['length'],width=r['width'],
        nup=r.get('nup'),cap=q['cap'],errors=errors,gates=gates,
        all_recorded_optimization_gates_passed=all(v is True for v in gates.values()),
        legacy_sweep_and_bond_gates_passed=legacy_passed,
        local_stationarity_evidence_scope="Sweep and bond criteria alone do not certify achieved local eigensolver accuracy or variational stationarity.",
        thermodynamic_length_convergence_certified=False,global_ground_sector_certified=False,
        common_topological_MES_certified=False,eligible_for_2D_phase_inference=False))
infinite=[]
for p in json.load(open('results/infinite_analysis.json'))['raw_measurements']:
    r=json.load(open(p['source_file']));audit=r.get('audit',{});iterations=audit.get('iterations',[])
    final=iterations[-1] if iterations else {}
    gates=dict(canonical_solver_residual=small(r.get('solver_residual'),thresholds['canonical_solver_residual']),
               consistent_measured_state=p['bound_comparison_state_consistency_checked'],
               all_local_eigenpairs_converged=final.get('all_local_eigenpairs_converged'))
    infinite.append(dict(source_file=p['source_file'],family=p['family'],width=p['width'],cap=p['cap'],
        allocated_chi=p['actual_chi'],gates=gates,
        all_recorded_optimization_gates_passed=all(v is True for v in gates.values()),
        circumference_convergence_certified=False,global_ground_density_certified=False,
        common_topological_MES_certified=False,physical_bulk_gap_certified=False,
        eligible_for_2D_phase_inference=False))
result=dict(thresholds=thresholds,threshold_source='Existing production controls plus the direct winding optimizer energy variance and known-sector ED gates. The demonstrated cold plateau requires separating legacy sweep/bond heuristics from stationarity; no tolerance has been relaxed.',
    null_gate_meaning='Required evidence is missing; unknown does not pass.',finite=finite,infinite=infinite,
    analysis_git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
    interpretation='Optimization accuracy and thermodynamic phase evidence are distinct. These gates describe recorded evidence; false phase-certification fields are not claims that the physical phase is absent.')
pathlib.Path('results/convergence_gate_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print('Optimization-gate passes:',sum(x['all_recorded_optimization_gates_passed'] for x in finite),'of',len(finite),'finite;',sum(x['all_recorded_optimization_gates_passed'] for x in infinite),'of',len(infinite),'infinite. No 2D phase certificate.')
