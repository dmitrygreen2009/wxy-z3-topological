"""Recover recorded iteration data without inventing missing launch provenance."""
import csv,json,pathlib,re
root=pathlib.Path('results')
finite=re.compile(r'After sweep (\d+) energy=([^ ]+)\s+maxlinkdim=(\d+) maxerr=([^ ]+) time=([^\s]+)')
infinite=re.compile(r'AUDITED VUMPS iteration=(\d+) chi=(\d+) residual=([^ ]+) seconds=([^\s]+)')
rows=[]
for path in sorted(root.glob('*.log')):
    for line_no,line in enumerate(path.read_text(errors='replace').splitlines(),1):
        m=finite.search(line)
        if m:
            sweep,e,chi,err,seconds=m.groups()
            rows.append(dict(log=str(path),line=line_no,solver='ITensorMPS DMRG',iteration=int(sweep),energy=float(e),bond_dimension=int(chi),truncation_error=float(err),canonical_residual=None,runtime_seconds=float(seconds)))
        m=infinite.search(line)
        if m:
            step,chi,err,seconds=m.groups()
            rows.append(dict(log=str(path),line=line_no,solver='ITensorInfiniteMPS VUMPS',iteration=int(step),energy=None,bond_dimension=int(chi),truncation_error=None,canonical_residual=float(err),runtime_seconds=float(seconds)))
with (root/'recorded_iterations.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
(root/'legacy_provenance_audit.json').write_text(json.dumps(dict(
    interpretation='Recorded log values, not reconstructed launch settings. Runtime excludes compilation, setup and measurements.',
    iteration_count=len(rows),
    original_parent_commit='80e59690c220e93b3241a9bf61d3eb1628da3956',
    first_reproducible_source_commit='5dfb1876638f8a38032c27f821174fc069849070',
    limitations=['Historical calculations preceded launch provenance. Original parent commit does not contain their uncommitted source.',
    'Seeds and settings must be read from preserved original scripts where known; missing values remain unknown.',
    'Original logs contain repeated stage counters; CSV identifies each value by log and line rather than asserting a unique run mapping.',
    'New version-3 jobs record source fingerprints and launch metadata in checkpoint sidecars.'],
    log_runtime_lower_bounds_seconds={str(p):sum(r['runtime_seconds'] for r in rows if r['log']==str(p)) for p in sorted(root.glob('*.log')) if any(r['log']==str(p) for r in rows)}),indent=2)+'\n')
print(f'Archived {len(rows)} recorded iterations; provenance gaps explicitly preserved.')
