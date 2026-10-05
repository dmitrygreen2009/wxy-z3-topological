"""Conservative SIGSTOP/SIGCONT control of an explicitly audited Julia queue.
Checks PID, launch time and command before signaling; never terminates a job.
"""
import argparse, datetime, json, os, pathlib, signal, subprocess, time


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
        for k in ['pid', 'ppid', 'started', 'command'])


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
    assert action == 'watch'
    print('Watching armchair winding batch; resume queue runs one heavy Julia job at a time.', flush=True)
    while True:
        current = processes()
        critical = data['critical_process']
        if same_process(critical, current.get(critical['pid'])) and 'Z' not in current[critical['pid']]['state']:
            time.sleep(30)
            continue
        complete = True
        for q in range(3):
            p = pathlib.Path(f'results/armchair_L2_w1_Nup8_winding{q}_qn_seed7254.json')
            if not p.exists() or not json.loads(p.read_text())['records'][-1]['converged_within_fixed_number_and_loop_sector']:
                complete = False
        if not complete:
            data['scheduler_status'] = 'Armchair batch exited before all three sector gates passed; queue remains paused'
            save(path, data)
            print(data['scheduler_status'], flush=True)
            return
        heavy = [p for p in current.values() if p['command'].startswith(data['julia_executable'])
                 and 'T' not in p['state'] and 'Z' not in p['state']]
        if heavy:
            time.sleep(30)
            continue
        resumed = False
        for job in data['resume_queue']:
            process = current.get(job['process']['pid'])
            if not same_process(job['process'], process) or 'Z' in process['state']:
                job['status'] = 'Process completed or no longer matches; no signal sent'
                continue
            if 'T' in process['state']:
                os.kill(process['pid'], signal.SIGCONT)
                job['status'] = 'SIGCONT sent; resumed without restart'
                job['resumed_utc'] = utc()
                print(utc(), 'Resumed', process['pid'], job['purpose'], flush=True)
                resumed = True
                break
        data['scheduler_status'] = 'One queued calculation resumed' if resumed else 'All queued calculations completed'
        save(path, data)
        if not resumed:
            return
        time.sleep(30)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['stop', 'watch'])
    parser.add_argument('manifest', type=pathlib.Path)
    args = parser.parse_args()
    control(args.manifest, args.action)
