"""Read ordinary staged runs and completed number-sector candidates without aliases."""
import glob,json,pathlib,re

def finite_groups(pattern=None):
    paths=glob.glob(pattern) if pattern is not None else glob.glob('results/*_L*_w*_chi*.json')+glob.glob('results/snapshots/*_L*_w*_chi*.json')+glob.glob('results/*_winding*_penalty_seed*.json')+glob.glob('results/*_winding*_qn_seed*.json')
    catalog=pathlib.Path('results/optimizer_fixture_catalog.json')
    fixtures=set(json.loads(catalog.read_text())['fixture_result_files']) if catalog.exists() else set()
    for path in sorted(set(paths)):
        r=json.loads(pathlib.Path(path).read_text());records=r.get('records',[])
        if r.get('validation_fixture',False) or path in fixtures:continue
        settings=r.get('audit',{}).get('solver_settings',{})
        if 'winding_charge' in settings:r['winding_charge']=settings['winding_charge']
        if not records:continue
        if 'family' in r and 'length' in r and all('maxlinkdim' in q for q in records):
            origin=r.get('snapshot_source_result',r.get('snapshot_audit',{}).get('source_result',path))
            r['initialization_branch']=re.sub(r'_chi\d+','',pathlib.Path(origin).stem).removesuffix('_refined')
            yield path,r,list(enumerate(records))
            continue
        # Number scans and PH-transformed checkpoints store their physical
        # metadata on each candidate instead of the outer envelope.
        for k,q in enumerate(records):
            if not {'family','length','width','physical_spins','nup','bond_dimension','cap','energy','entropy'}<=q.keys():continue
            meta=dict(q);meta['spins']=q['physical_spins']
            meta['initialization_branch']=re.sub(r'_chi\d+','',pathlib.Path(path).stem)+f'_Nup{q["nup"]}_seed{q.get("seed","unknown")}'
            stage=dict(q);stage['maxlinkdim']=q['bond_dimension']
            stage['sweep_max_truncation_errors']=q.get('truncation_errors') or [None]
            yield path,meta,[(k,stage)]
