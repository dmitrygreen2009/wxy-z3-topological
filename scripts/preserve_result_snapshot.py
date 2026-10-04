"""Freeze a completed result and its local checkpoint before another run reuses an alias."""
import argparse, datetime, hashlib, json, pathlib, shutil, subprocess
p = argparse.ArgumentParser()
p.add_argument('result', type=pathlib.Path)
p.add_argument('label')
a = p.parse_args()
if not a.label.replace('_', '').replace('-', '').isalnum():
    raise ValueError('Use a plain snapshot label')
meta = json.loads(a.result.read_text())
payload = pathlib.Path(meta['checkpoint_file'])
def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()
source_hash = digest(payload)
if 'checkpoint_sha256' in meta and source_hash != meta['checkpoint_sha256']:
    raise ValueError('Result checksum does not match checkpoint')
stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
copy = payload.with_name(payload.stem + '_' + a.label + '_' + stamp + payload.suffix)
shutil.copy2(payload, copy)
if digest(copy) != source_hash:
    raise ValueError('Snapshot checksum mismatch; original preserved')
sidecar = payload.with_suffix('.json')
if sidecar.exists():
    checkpoint_meta = json.loads(sidecar.read_text())
    checkpoint_meta['checkpoint_file'] = str(copy)
    checkpoint_meta['checkpoint_sha256'] = source_hash
    copy.with_suffix('.json').write_text(json.dumps(checkpoint_meta, indent=2) + '\n')
meta['checkpoint_file'] = str(copy)
meta['checkpoint_sha256'] = source_hash
meta['checkpoint_bytes'] = copy.stat().st_size
meta['snapshot_audit'] = {'source_result': str(a.result), 'source_result_sha256': digest(a.result),
    'source_checkpoint': str(payload), 'created_utc': stamp,
    'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'operation': 'Exact byte copy; no optimization or numerical observable recomputation'}
out = pathlib.Path('results/snapshots') / (a.label + '_' + stamp + '.json')
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(meta, indent=2) + '\n')
print(out)
