"""Conservative SIGSTOP/SIGCONT control of an explicitly audited Julia queue.
Checks PID, launch time and command before signaling; never terminates a job.
"""
import argparse, datetime, fcntl, hashlib, json, math, os, pathlib, re, signal, subprocess, time


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def processes():
    output = subprocess.check_output(
        ['ps', '-axo', 'pid=,ppid=,stat=,lstart=,command='], text=True)
    result = {}
    for line in output.splitlines():
        parts = line.split(None, 8)
        if len(parts) == 9:
            pid, ppid, state = parts[:3]
            result[int(pid)] = dict(pid=int(pid), ppid=int(ppid), state=state,
                started=' '.join(parts[3:8]), command=parts[8])
    return result


def same_process(saved, current):
    return current is not None and all(saved[k] == current[k]
        for k in ['pid', 'started', 'command'])


def save(path, data):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    os.replace(temporary, path)


def control(path, action):
    data = json.loads(path.read_text())
    if action == 'stop':
        for job in data['resume_queue']:
            current = processes().get(job['process']['pid'])
            if not same_process(job['process'], current):
                job['status'] = 'Exited or identity changed; no signal sent'
            else:
                assert pathlib.Path(job['checkpoint_metadata']).is_file()
                assert pathlib.Path(job['checkpoint_payload']).is_file()
                os.kill(current['pid'], signal.SIGSTOP)
                job['status'] = 'SIGSTOP sent; in-memory state and parent retained'
                job['paused_utc'] = utc()
            save(path, data)
        return
    watch(path)


def announce(path, data, state, reason):
    changed = data.get('queue_state') != state or data.get('scheduler_status') != reason
    data.update(queue_state=state, scheduler_status=reason, scheduler_checked_utc=utc())
    if changed:
        message = f"{utc()} {state}: {reason}"
        print(message, flush=True)
        if state in ('BLOCKED', 'COMPLETE'):
            marker = path.parent / f'QUEUE_{state}.txt'
            marker.write_text(message + '\n')
            try:
                script = 'display notification ' + json.dumps(reason[:240]) + ' with title "WXY queue ' + state + '"'
                subprocess.run(['osascript', '-e', script], timeout=5, capture_output=True)
            except (OSError, subprocess.TimeoutExpired):
                pass
    if state != 'BLOCKED':
        (path.parent / 'QUEUE_BLOCKED.txt').unlink(missing_ok=True)
    if state != 'COMPLETE':
        (path.parent / 'QUEUE_COMPLETE.txt').unlink(missing_ok=True)
    save(path, data)


def checkpoint_record(job):
    # Record artifacts only; solver exit is not evidence of physical convergence.
    result = {'source_checkpoint': job.get('checkpoint_payload'), 'log': job.get('restart_log', job.get('log'))}
    plan = job.get('restart_plan')
    if plan:
        stem = plan['snapshot'].split('_latest')[0]
        files = list(pathlib.Path(stem).parent.glob(re.sub(r'_chi\d+_', '_chi*_', pathlib.Path(stem).name) + '*.json'))
        files = [p for p in files if 'resource_restart_snapshot' not in p.name]
        if files:
            latest = max(files, key=lambda p: p.stat().st_mtime)
            metadata = json.loads(latest.read_text())
            result.update(latest_metadata=str(latest), latest_checkpoint=metadata.get('checkpoint_file'),
                          final_cap=metadata.get('cap'), iteration=metadata.get('iteration'))
    else:
        files = list(pathlib.Path('results/checkpoints').glob('*.json'))
        if files:
            latest = max(files, key=lambda p: p.stat().st_mtime)
            result['latest_checkpoint_metadata_at_exit'] = str(latest)
            result['checkpoint_scope'] = 'Most recently written metadata; not a certificate of this batch completion'
    return result


def inspect_completion(job):
    """Check the output contract; scheduled completion is distinct from convergence."""
    expected = job.get('expected_result')
    if expected:
        result_path = pathlib.Path(expected['path'])
        if not result_path.exists():
            return False, {'error': 'Expected result is missing: ' + str(result_path)}
        metadata = json.loads(result_path.read_text())
        for key, value in expected.get('fields', {}).items():
            if metadata.get(key) != value:
                return False, {'error': f'Expected {key}={value}, got {metadata.get(key)} in {result_path}'}
        if 'final_phase' in expected and metadata.get('phase') != expected['final_phase']:
            return False, {'error': 'Incomplete final phase: ' + str(metadata.get('phase'))}
        records = metadata.get('records', [])
        if expected.get('passes'):
            tail = records[-expected['passes']:]
            if len(tail) != expected['passes'] or any(r.get('cap') != expected['fields']['cap'] or len(r.get('sweep_energies', [])) != expected['sweeps_per_pass'] for r in tail):
                return False, {'error': 'Required scheduled passes/sweeps missing'}
    elif job.get('restart_plan'):
        artifacts = checkpoint_record(job)
        if not artifacts.get('latest_metadata'):
            return False, {'error': 'No final checkpoint metadata for restart'}
        metadata = json.loads(pathlib.Path(artifacts['latest_metadata']).read_text())
        plan = job['restart_plan']
        if metadata.get('cap') != plan['target_cap']:
            return False, {'error': 'Target bond-dimension stage not completed'}
        # An atomic stage-complete manifest is written only after measurement.
        stem = plan['snapshot'].split('_latest')[0]
        # Stage-complete filenames also include Ly; match metadata, not implicit paths.
        source = json.loads(pathlib.Path(plan['snapshot'].replace('.jls', '.json')).read_text())
        complete = []
        pattern = f"infinite_{source['family']}_w{source['width']}_*stage*_complete.json"
        for p in pathlib.Path('results/checkpoints').glob(pattern):
            m = json.loads(p.read_text())
            fields = ('family', 'width', 'cell_spins', 'infinite_ordering', 'measurement_tag')
            if m.get('cap') == plan['target_cap'] and all(m.get(k) == source.get(k) for k in fields):
                complete.append((p, m))
        if not complete:
            return False, {'error': 'Target stage completion manifest missing'}
        result_path, metadata = max(complete, key=lambda pair: pair[0].stat().st_mtime)
    else:
        return False, {'error': 'No explicit expected-output contract; completion is ambiguous'}
    checkpoint = pathlib.Path(metadata['checkpoint_file'])
    sidecar = pathlib.Path(str(checkpoint).replace('.jls', '.json'))
    cp = json.loads(sidecar.read_text())
    if not checkpoint.exists() or checkpoint.stat().st_size != cp['checkpoint_bytes']:
        return False, {'error': 'Checkpoint missing or size mismatch: ' + str(checkpoint)}
    with checkpoint.open('rb') as f:
        digest = hashlib.file_digest(f, 'sha256').hexdigest()
    if digest != cp['checkpoint_sha256'] or digest != metadata['checkpoint_sha256']:
        return False, {'error': 'Checkpoint checksum mismatch: ' + str(checkpoint)}
    if job.get('restart_launch_utc'):
        written = datetime.datetime.fromisoformat(cp['checkpoint_written_utc'].replace('Z', '+00:00'))
        if written.tzinfo is None:
            written = written.replace(tzinfo=datetime.timezone.utc)
        launched = datetime.datetime.fromisoformat(job['restart_launch_utc'])
        if written < launched:
            return False, {'error': 'Completion checkpoint predates this launch'}
    scalar = metadata.get('solver_residual')
    energy = metadata.get('energy_per_vertex', metadata.get('energy'))
    if metadata.get('records'):
        energy = metadata['records'][-1]['energy']
    if energy is not None and not math.isfinite(energy):
        return False, {'error': 'Nonfinite final energy'}
    if job.get('requires_solver_convergence') and (scalar is None or scalar > metadata.get('tol', 1e-7)):
        return False, {'error': 'Required original solver convergence criterion failed'}
    return True, {'result': str(result_path), 'checkpoint': str(checkpoint), 'sha256': digest,
                  'scheduled_completion_validated': True, 'energy': energy,
                  'solver_residual': scalar, 'phase_convergence_certified': False,
                  'interpretation': 'Scheduled computation finished; existing scientific convergence gates remain unchanged'}


def record_completion(job):
    valid, evidence = inspect_completion(job)
    job['completion_inspection'] = evidence
    if valid:
        job.update(completion_validated=True, execution_finished=True, completed_utc=utc())
    return valid, evidence


def tick(path, children):
    data = json.loads(path.read_text())
    current = processes()
    if data.get('queue_state') == 'BLOCKED' and not data.get('blocked_job_pid'):
        return
    # Collect genuine exit status for children owned by this watcher.
    for job in data['resume_queue']:
        key = str(job['process']['pid'])
        child = children.get(key)
        if child is not None and child.poll() is not None and not job.get('execution_finished'):
            job.update(execution_finished=True, exit_code=child.returncode, exited_utc=utc(),
                       exit_evidence='waitpid via Popen.poll', final_artifacts=checkpoint_record(job))
            print(utc(), 'EXIT', job['purpose'], 'status', child.returncode, json.dumps(job['final_artifacts']), flush=True)
            valid, evidence = record_completion(job)
            if not valid:
                data['blocked_job_pid'] = job['process']['pid']
                announce(path, data, 'BLOCKED', job['purpose'] + ': ' + evidence['error'])
                return
    heavy = [p for p in current.values() if p['command'].startswith(data['julia_executable'])
             and 'T' not in p['state'] and 'Z' not in p['state']]
    if heavy:
        announce(path, data, 'RUNNING', 'Active Julia PID(s): ' + ','.join(str(p['pid']) for p in heavy))
        return
    if data.get('intentional_blocker'):
        announce(path, data, 'BLOCKED', str(data['intentional_blocker']))
        return
    critical = data['critical_process']
    if same_process(critical, current.get(critical['pid'])) and 'Z' not in current[critical['pid']]['state']:
        announce(path, data, 'BLOCKED', 'Critical process is stopped; requires explicit resume authorization: ' + str(critical['pid']))
        return
    for q in range(3):
        p = pathlib.Path(f'results/armchair_L2_w1_Nup8_winding{q}_qn_seed7254.json')
        if not p.exists() or not json.loads(p.read_text())['records'][-1]['converged_within_fixed_number_and_loop_sector']:
            announce(path, data, 'BLOCKED', f'Existing armchair winding gate q={q} failed/missing; inspect results before advancing')
            return
    for job in data['resume_queue']:
        if job.get('execution_finished'):
            if data.get('blocked_job_pid') == job['process']['pid']:
                valid, evidence = record_completion(job)
                if not valid:
                    announce(path, data, 'BLOCKED', job['purpose'] + ': ' + evidence['error'])
                    return
                data.pop('blocked_job_pid', None)
            continue
        saved = job.get('restart_process') if job.get('restart_launched') else job['process']
        process = current.get(saved['pid']) if saved else current.get(job.get('restart_pid'))
        if job.get('restart_launched') and saved is None and process is not None:
            # Adopt a pre-fix live restart only if its checkpoint identity agrees.
            if job['restart_plan']['snapshot'] not in process['command']:
                announce(path, data, 'BLOCKED', 'Restart PID identity changed; refusing duplicate: ' + job['purpose'])
                return
            saved = process
            job['restart_process'] = process
        alive = saved is not None and same_process(saved, process) and 'Z' not in process['state']
        if alive:
            if 'T' in process['state']:
                announce(path, data, 'TRANSITIONING', 'Resuming preserved in-memory PID ' + str(process['pid']))
                os.kill(process['pid'], signal.SIGCONT)
                job.update(resumed_utc=utc(), status='SIGCONT: existing in-memory solver state preserved',
                           protected_resume_authorization='Explicit user request 2026-10-07')
                announce(path, data, 'RUNNING', job['purpose'] + ' PID ' + str(process['pid']))
                return
        if job.get('restart_launched') or not job.get('restart_ready'):
            job.update(execution_finished=True, exit_code=None, exited_utc=utc(),
                       exit_evidence='Pre-existing non-child process disappeared; exit status unavailable, success NOT inferred',
                       final_artifacts=checkpoint_record(job))
            print(utc(), 'EXIT STATUS UNKNOWN', job['purpose'], json.dumps(job['final_artifacts']), flush=True)
            if job.get('protected') or job.get('restart_launched'):
                valid, evidence = record_completion(job)
                if not valid:
                    data['blocked_job_pid'] = job['process']['pid']
                    announce(path, data, 'BLOCKED', job['purpose'] + ': ' + evidence['error'])
                    return
                print(utc(), 'OUTPUT VALIDATED', job['purpose'], json.dumps(evidence), flush=True)
            continue
        if not job.get('restart_released') or not job.get('restart_validation', {}).get('validated'):
            announce(path, data, 'BLOCKED', 'No validated/released restart path for ' + job['purpose'])
            return
        plan = job['restart_plan']
        snapshot = pathlib.Path(plan['snapshot'])
        with snapshot.open('rb') as source:
            if hashlib.file_digest(source, 'sha256').hexdigest() != plan['snapshot_sha256']:
                announce(path, data, 'BLOCKED', 'Checkpoint hash mismatch; restore validated payload: ' + str(snapshot))
                return
        # Single watcher lock and fresh active-process check prohibit duplicates.
        if any(p['command'].startswith(data['julia_executable']) and 'T' not in p['state'] and 'Z' not in p['state'] for p in processes().values()):
            announce(path, data, 'TRANSITIONING', 'A Julia job became active; defer launch to next check')
            return
        command = [data['julia_executable'], '--project=.', '-e',
            'include("scripts/resume_checkpoint.jl");resume_checkpoint(ARGS[1],parse(Int,ARGS[2]);preserve_stage_budget=true)',
            plan['snapshot'], str(plan['target_cap'])]
        environment = os.environ.copy()
        environment.update(JULIA_DEPOT_PATH='/private/tmp/wxy-julia-depot', OPENBLAS_NUM_THREADS='1')
        announce(path, data, 'TRANSITIONING', 'Launching validated checkpoint: ' + job['purpose'])
        with open(job['restart_log'], 'ab') as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=environment, start_new_session=True)
        children[str(job['process']['pid'])] = child
        job.update(restart_launched=True, restart_pid=child.pid, restart_launch_utc=utc(), restart_command=command,
                   restart_process=processes().get(child.pid))
        announce(path, data, 'RUNNING', job['purpose'] + ' PID ' + str(child.pid))
        return
    unverified = [j['purpose'] for j in data['resume_queue'] if not (j.get('completion_validated') or j.get('scientific_completion_verified'))]
    if unverified:
        announce(path, data, 'BLOCKED', 'All queued processes exited, but completion/convergence requires result validation: ' + '; '.join(unverified))
    else:
        announce(path, data, 'COMPLETE', json.dumps([{'job': j['purpose'], 'artifacts': j.get('final_artifacts')} for j in data['resume_queue']]))


def watch(path):
    # OS releases this simple lock if the watcher dies. No duplicate watchers.
    with open(str(path) + '.lock', 'w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        children = {}
        while True:
            try:
                tick(path, children)
            except Exception as error:
                data = json.loads(path.read_text())
                announce(path, data, 'BLOCKED', f'Watcher recovery failed: {type(error).__name__}: {error}; correct queue/checkpoint access')
            time.sleep(30)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['stop', 'watch'])
    parser.add_argument('manifest', type=pathlib.Path)
    args = parser.parse_args()
    control(args.manifest, args.action)
