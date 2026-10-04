"""Exact-spin number-sector lower bounds; reuses saved validation results."""
import datetime,hashlib,json,math,pathlib,subprocess
from fractions import Fraction
started=datetime.datetime.now(datetime.timezone.utc)
root=pathlib.Path(__file__).resolve().parents[1]
flat=json.loads((root/'results/flatband_audit.json').read_text())
assert all(p['maximum_error']<1e-12 for p in flat['points'])
ed=json.loads((root/'results/ed.json').read_text());checks=[]
for name,r in ed.items():
    legs=r['legs_zero_based'];gauges=max(e for row in legs for e in row)+1
    spins=3*len(legs)+gauges
    degree=max(sum(e in row for row in legs) for e in range(gauges))
    lower=-math.sqrt(degree)*min(r['nup'],spins-r['nup'])
    assert min(r['energies'])>=lower-1e-10
    checks.append(dict(benchmark=name,spins=spins,nup=r['nup'],maximum_edge_degree=degree,
        exact_sector_lower_bound=lower,saved_ground_energy=min(r['energies'])))
# Valid in the actual 2D geometry: disjoint A stars, fixed-number B matter.
# Round the validated cluster energy upward to retain a conservative trial.
trial_per_primitive_cell=-2.4030921
minimum_density=-trial_per_primitive_cell/(9*math.sqrt(2))
grid=[Fraction(k,27) for k in range(28)]
remaining=[str(q) for q in grid if minimum_density<=q<=1-minimum_density]
assert remaining==[str(Fraction(k,27)) for k in range(6,22)]
r=dict(derivation='The microscopic one-particle hopping matrix h has minimum eigenvalue -sqrt(d_max), since W is unitary and C dagger C is the diagonal physical-edge incidence count. Factoring h+sqrt(d_max) I as L dagger L gives H+sqrt(d_max) N_up=sum_a B_a dagger B_a >=0 with B_a=sum_j L_aj S_j minus. Applying the same argument to h transpose and raising operators gives H>=-sqrt(d_max)(N_spins-N_up). Therefore E(N)>=-sqrt(d_max) min(N,N_spins-N). Hard-core spin algebra is retained; no free-particle many-body replacement is made.',
    benchmark_checks=checks,many_body_diagonalizations_repeated=False,
    bulk_trial_energy_per_primitive_cell=trial_per_primitive_cell,
    bulk_density_window_from_exact_lower_bound_and_trial=[minimum_density,1-minimum_density],
    ordinary_D_Z3_density_grid_remaining_after_energy_bound=remaining,
    interpretation='This excludes extreme densities as global ground candidates using an exact operator inequality and a valid 2D trial. It does not select the preferred density or identify a topological phase. A stronger trial on a narrow cylinder must not be used as a 2D bulk energy bound without constructing a valid 2D extension.',
    source_validation_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'results/flatband_audit.json',root/'results/ed.json']},
    git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
    analyzed_utc=started.isoformat())
(root/'results/number_bound_audit.json').write_text(json.dumps(r,indent=2)+'\n')
print('Exact number bound checked against preserved benchmarks; ordinary D(Z3) density candidates remaining:',len(remaining))
