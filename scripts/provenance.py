import hashlib,pathlib,platform,subprocess,datetime

def provenance(seed,solver,settings,initialization,conserved):
    try:
        commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
        dirty=bool(subprocess.check_output(['git','status','--porcelain','--','src','scripts','test','Project.toml','Manifest.toml','requirements.txt','requirements.lock.txt'],text=True).strip())
    except Exception:commit='unknown';dirty=True
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for folder in ['src','scripts','test'] for p in pathlib.Path(folder).iterdir() if p.is_file()}
    return dict(random_seed=seed,solver=solver,settings=settings,initialization=initialization,
        conserved_quantum_numbers=conserved,git_commit=commit,code_worktree_dirty=dirty,
        source_sha256=hashes,python_version=platform.python_version(),started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
