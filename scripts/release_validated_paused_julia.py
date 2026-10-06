"""Release only independently validated, paused checkpoint-restart jobs.
Uses SIGTERM then SIGCONT to deliver termination to a stopped process.
Never SIGKILLs or signals the active priority calculation.
"""
import hashlib, json, os, pathlib, signal, subprocess, time
from julia_resource_queue import processes, same_process, save, utc

path = pathlib.Path('results/julia_resource_pause_audit.json')
verification = json.loads(pathlib.Path('results/julia_checkpoint_restart_validation.json').read_text())
verified = {r['pid']:r for r in verification['records'] if r['validated']}
plans = json.loads(pathlib.Path('results/julia_checkpoint_restart_plans.json').read_text())
audit_path = pathlib.Path('results/julia_checkpoint_ram_release_audit.json')
record = json.loads(audit_path.read_text()) if audit_path.exists() else {'started_utc':utc(), 'records':[], 'protected_pid':plans['protected_finite_job']['pid']}
record['latest_release_attempt_utc'] = utc()

def memory():
    result = {}
    for name, command in [('vm_stat',['vm_stat']), ('swap_usage',['sysctl','-n','vm.swapusage']), ('memory_pressure',['memory_pressure','-Q'])]:
        p = subprocess.run(command, capture_output=True, text=True)
        result[name] = {'returncode':p.returncode, 'stdout':p.stdout, 'stderr':p.stderr}
    return result

if 'memory_before_release' not in record:
    record['memory_before_release'] = memory()
for plan in plans['jobs']:
    pid = plan['pid']; data = json.loads(path.read_text())
    job = next(j for j in data['resume_queue'] if j['process']['pid']==pid)
    assert pid in verified and job['restart_ready'] and pid!=data['active_priority_pid_at_conversion']
    if job.get('restart_released'):
        continue
    pending = any(r['pid']==pid and not r['termination_confirmed'] for r in record['records'])
    current = processes().get(pid)
    already_exited = not same_process(plan['process'], current)
    assert already_exited and pending or same_process(plan['process'],current) and 'T' in current['state'],pid
    with open(plan['snapshot'],'rb') as f:
        assert hashlib.file_digest(f,'sha256').hexdigest()==plan['snapshot_sha256']
    if not already_exited:
        if not pending:
            os.kill(pid,signal.SIGTERM)
        try: os.kill(pid,signal.SIGCONT)
        except ProcessLookupError: pass
    exited = False
    for _ in range(480):
        current = processes().get(pid)
        if not same_process(plan['process'],current) or 'Z' in current['state']:
            exited = True;break
        time.sleep(.25)
    if not exited:
        os.kill(pid,signal.SIGSTOP)
    job['restart_released'] = exited
    job['restart_release_utc'] = utc()
    job['restart_status'] = 'Original paused process exited after SIGTERM; immutable restart snapshot retained' if exited else 'Termination not confirmed; process stopped again, no automatic replacement launch'
    job['status'] = job['restart_status']
    save(path,data)
    record['records'] = [r for r in record['records'] if r['pid']!=pid]
    record['records'].append(dict(pid=pid,termination_confirmed=exited,snapshot=plan['snapshot'],checkpoint_sha256=plan['snapshot_sha256'],unfinished_work=plan['unfinished_work'],completed_unsaved_iterations=plan['completed_unsaved_iterations']))
    save(pathlib.Path('results/julia_checkpoint_ram_release_audit.json'),record)
    print(pid,job['restart_status'],flush=True)
record['memory_immediately_after_release'] = memory()
record['finished_utc'] = utc()
save(pathlib.Path('results/julia_checkpoint_ram_release_audit.json'),record)
