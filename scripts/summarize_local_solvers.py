"""Snapshot complete local-solver evidence without changing a live calculation."""
import argparse,datetime,hashlib,json,pathlib,subprocess
parser=argparse.ArgumentParser()
parser.add_argument('source',type=pathlib.Path)
parser.add_argument('output',type=pathlib.Path)
parser.add_argument('--through',type=int)
args=parser.parse_args()
data=args.source.read_bytes();prefix=[];records=[]
for line in data.splitlines(keepends=True):
    if not line.endswith(b'\n'):
        break  # A writer may still be appending the last line.
    r=json.loads(line)
    if args.through is not None and r['iteration']>args.through:
        break
    solves=r['solves'];prefix.append(line)
    records.append(dict(iteration=r['iteration'],git_commit=r.get('git_commit'),
        local_solve_count=len(solves),all_local_criteria_passed=bool(solves) and all(s['local_criterion_passed'] for s in solves),
        maximum_true_residual=max((s['true_eigenpair_residual'] for s in solves),default=None)))
snapshot=args.output.with_suffix('.jsonl');snapshot.parent.mkdir(parents=True,exist_ok=True)
snapshot.write_bytes(b''.join(prefix))
result=dict(audited_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_jsonl=str(args.source),
    preserved_source_prefix=str(snapshot),preserved_source_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),
    analysis_git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),records=records,
    interpretation='Complete logged local-solve evidence only. Local accuracy does not certify canonical self-consistency, global ground energy, or thermodynamic convergence. Source calculation is unchanged.')
args.output.write_text(json.dumps(result,indent=2)+'\n')
print(len(records),'completed iterations; all local criteria passed:',bool(records) and all(r['all_local_criteria_passed'] for r in records))
