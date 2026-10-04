"""Independent exact N_up=1 spectrum; no many-body phase inference."""
import json,pathlib,time
import numpy as np
from ed import hamiltonian
from provenance import provenance
paths=[pathlib.Path('geometry')/name for name in ['single_star.json','shared_edge_pair.json','torus_four_star.json']]
paths += [pathlib.Path(f'geometry/{family}_L2_w{w}_star.json') for family in ['zigzag','armchair'] for w in [1,2,3]]
results=[]
for path in paths:
    started=time.time();g=json.loads(path.read_text());legs=[[e-1 for e in es] for es in g['endpoint_leg_table']]
    h,_,herm=hamiltonian(legs,1)
    actual=np.linalg.eigvalsh(h.toarray())
    multiplicities=np.array([len(e['endpoints']) for e in g['shared_gauge_edges']])
    nm=len(g['matter_spins']);ng=len(multiplicities)
    expected=np.sort(np.r_[-np.sqrt(multiplicities),np.zeros(nm-ng),np.sqrt(multiplicities)])
    error=float(max(abs(actual-expected)))
    assert herm<1e-13 and error<1e-12
    results.append(dict(geometry=str(path),spins=len(actual),nup=1,energies=actual.tolist(),expected_energies=expected.tolist(),maximum_error=error,hermiticity_error=herm,runtime_seconds=time.time()-started,
        audit=provenance(None,'Independent dense exact one-particle spin-sector ED',dict(tolerance=1e-12),'Full exact N_up=1 matrix',['N_up=1'])))
pathlib.Path('results/flatband_audit.json').write_text(json.dumps(dict(points=results,
    derivation='W dagger W=I makes the matter-to-gauge hopping matrix satisfy A dagger A=diag(edge endpoint counts). Exact N_up=1 energies are paired +/-sqrt(count), plus Nm-Ng zeros.',
    interpretation='Validates microscopic shared-edge incidence on both cylinder families. Flat single-particle bands do not determine the interacting many-body ground-state phase.'),indent=2)+'\n')
print(f'Passed exact one-particle spectra for {len(results)} geometries; max error {max(x["maximum_error"] for x in results):.3g}.')
