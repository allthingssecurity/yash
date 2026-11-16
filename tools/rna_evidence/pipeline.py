"""Plan, run, collect or replay reusable pre-selection RNA evidence screening."""
from pathlib import Path
import argparse,csv,gzip,hashlib,io,json,re,sys
ROOT=Path(__file__).resolve().parents[2]
IMAGE='im-hX1h6HknDMvasm9Z6q530H'
VOLUME='research-lens-workspace'
PARENT='/data/results/TN25-112103/matched-normal-20261107'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def worker_bundle():
    return [ROOT/'tools/rna_evidence/pool_worker.py',ROOT/'tools/rna_evidence/worker.py',ROOT/'designs/design-2/pipeline/normal_evidence.py',ROOT/'designs/design-2/pipeline/self_search.py']


def zipped(path,value):
    body=(json.dumps(value,separators=(',',':'),sort_keys=True)+'\n').encode()
    buf=io.BytesIO()
    with gzip.GzipFile(filename='',mode='wb',fileobj=buf,mtime=0) as f:f.write(body)
    return buf.getvalue()


def prepare(candidates):
    source=ROOT/'data/shared-rna/splicing';normal=source/'normal-screen'
    with gzip.open(normal/'all-windows.tsv.gz','rt') as f:baseline={r['window_id']:r for r in csv.DictReader(f,delimiter='\t')}
    if candidates:
        path=(ROOT/candidates).resolve()
        if path.suffix=='.json':requested={r['window_id'] for r in json.loads(path.read_text()) if r.get('source')=='SPLICE'}
        else:
            with (gzip.open(path,'rt') if str(path).endswith('.gz') else path.open()) as f:requested={r['window_id'] for r in csv.DictReader(f,delimiter='\t')}
        assert requested and requested<=set(baseline),'Unknown RNA window IDs'
    else:requested={wid for wid,r in baseline.items() if r['primary_window_eligible']=='True'}
    assert all(baseline[wid]['primary_window_eligible']=='True' for wid in requested),'This stage accepts reconstructed windows passing the recorded earlier QC/normal/strand/HLA gates'
    paths={'reconstruction':source/'reconstruction/patient-windows.jsonl.gz','RNA_BAM_verification':source/'read-verification/verification.json','normal_screen':normal/'all-windows.tsv.gz','normal_screen_verification':normal/'verification.json','normal_selection_policy':normal/'selection-policy.json','affinities':source/'scoring/affinities.tsv.gz','native_reference_matches':source/'reference-check/matches.json.gz'}
    verification=json.loads(paths['RNA_BAM_verification'].read_text());assert verification['status']=='PASS' and verification['reconstruction_sha256']==sha(paths['reconstruction'])
    assert json.loads(paths['normal_screen_verification'].read_text())['status']=='PASS'
    windows=[]
    with gzip.open(paths['reconstruction'],'rt') as f:
        for line in f:
            w=json.loads(line)
            if w['window_id'] not in requested:continue
            r=baseline[w['window_id']];assert r['sequence']==w['sequence'] and r['normal_status']=='PASS' and r['strand_status']=='PASS'
            w['region_id']=r['region_id'];w['baseline']={k:r[k] for k in ['native_affinity_sum','weighted_score_left','seed_peptides','best_target_affinity_nm','normal_status','normal_detection','strand_status','strand_agreeing_fragments','strand_opposing_fragments']}
            windows.append(w)
    assert len(windows)==len(requested)
    actual_peptides={w['sequence'][i:i+k] for w in windows for k in [8,9,10,11] for i in range(26-k) if i<w['target_end_in_window'] and i+k>w['target_start_in_window']}
    affinities={p:{} for p in actual_peptides}
    with gzip.open(paths['affinities'],'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['peptide'] in affinities:
                assert r['allele'] not in affinities[r['peptide']]
                affinities[r['peptide']][r['allele']]=float(r['affinity_nm'])
    alleles=json.loads((source/'policy.json').read_text())['hla_alleles']
    assert len(alleles)==6 and all(set(v)==set(alleles) for v in affinities.values())
    with gzip.open(paths['native_reference_matches'],'rt') as f:native=json.load(f)
    policy_path=ROOT/'tools/rna_evidence/pool-policy.json';policy=json.loads(policy_path.read_text())
    u=ROOT/'designs/design-3/inputs/sequence-review/uniprot-reference-check.json';uniprot=json.loads(u.read_text())['reference_download']
    paths.update(screening_policy=policy_path,uniprot_source=u)
    return dict(schema_version=1,scope='RNA windows before additional antigen selection; no automatic construct change',windows=windows,policy=policy,alleles=alleles,
      affinities=affinities,original_native_reference_matches={p:bool(native[p]['occurs']) for p in sorted(actual_peptides)},
      reference_ensembl_sha256='f8df3d7cfc85718e9161579bd9168a3598547ee50022b5e1c272ed66293f0505',reference_uniprot=uniprot,
      source_hashes={k:dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for k,p in paths.items()})


def verify_outputs(out):
    receipt=json.loads((out/'complete.json').read_text());assert receipt['status']=='COMPLETED'
    assert receipt['request_sha256']==sha(out/'request.json.gz')
    bundle=worker_bundle()
    assert receipt['worker_sha256']==sha(bundle[0]), 'Preserve the measured worker source'
    state=json.loads((out/'state.json').read_text())
    digest=hashlib.sha256(b''.join(p.name.encode()+p.read_bytes() for p in bundle)).hexdigest()
    assert state['worker_bundle_sha256']==digest and state['request_sha256']==receipt['request_sha256']
    for name,r in receipt['outputs'].items():
        p=out/name;assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],name
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates',type=Path,help='Optional JSON/TSV window list; default is every previously eligible ND/strand-PASS RNA window')
    p.add_argument('--output',type=Path,default=Path('work/rna-pool-screen-20261112'))
    p.add_argument('--run-id',default='rna-pool-screen-20261112')
    p.add_argument('--check-selection',type=Path,help='With --verify, require proposed SPLICE-source regions to be in the primary checked export before assembly')
    mode=p.add_mutually_exclusive_group();mode.add_argument('--execute',action='store_true');mode.add_argument('--collect',action='store_true');mode.add_argument('--verify',action='store_true')
    a=p.parse_args();assert re.fullmatch('[a-z0-9-]+',a.run_id)
    if a.check_selection and not a.verify:p.error('--check-selection requires --verify')
    out=(ROOT/a.output).resolve()
    if a.verify:
        verify_outputs(out)
        from pool_report import replay_and_report
        print(json.dumps(replay_and_report(out),indent=2))
        if a.check_selection:
            from pool_report import require_selection
            print(json.dumps(require_selection(out,(ROOT/a.check_selection).resolve()),indent=2))
            from normal_peptides import require_selection as require_nonself
            print(json.dumps(require_nonself((ROOT/a.check_selection).resolve(), screen=out),indent=2))
        return
    request=prepare(a.candidates);body=zipped(None,request)
    frozen_request=out/'request.json.gz'
    if frozen_request.exists():
        with gzip.open(frozen_request,'rt') as f:recorded=json.load(f)
        assert recorded==request,'Different logical inputs require a new output directory'
        body=frozen_request.read_bytes()
    request_sha=hashlib.sha256(body).hexdigest()
    bundle=worker_bundle();worker=bundle[0]
    digest=hashlib.sha256(b''.join(p.name.encode()+p.read_bytes() for p in bundle)).hexdigest()
    remote=f'/data/results/TN25-112103/{a.run_id}/{request_sha[:12]}-{digest[:12]}'
    print(json.dumps(dict(windows=len(request['windows']),events=len({w['event_id'] for w in request['windows']}),genes=len({w['gene_id'] for w in request['windows']}),prediction_pairs=len(request['affinities'])*len(request['alleles']),new_model_inference=False,remote=remote,execute=a.execute,collect=a.collect),indent=2),flush=True)
    if not a.execute and not a.collect:return
    import modal
    out.mkdir(parents=True,exist_ok=True);req=out/'request.json.gz';statefile=out/'state.json'
    if req.exists():assert req.read_bytes()==body,'Different inputs require a new output directory'
    else:req.write_bytes(body)
    state=json.loads(statefile.read_text()) if statefile.exists() else dict(remote=remote,request_sha256=request_sha,worker_bundle_sha256=digest)
    assert state['remote']==remote and state['worker_bundle_sha256']==digest
    volume=modal.Volume.from_name(VOLUME)
    if 'sandbox_id' not in state:
        if a.collect:raise RuntimeError('Run --execute first')
        uploads=[(req,remote+'/request.json.gz')]+[(p,remote+'/code/'+p.name) for p in bundle]
        with volume.batch_upload() as upload:
            for path,key in uploads:upload.put_file(str(path),key.removeprefix('/data/'))
        sandbox=modal.Sandbox.create('python','-u',remote+'/code/pool_worker.py','--root',PARENT,'--request',remote+'/request.json.gz','--output',remote+'/results',
          app=modal.App.lookup('research-yash-matched-normal'),image=modal.Image.from_id(IMAGE),cpu=4,memory=24576,timeout=10800,volumes={'/data':volume})
        state.update(sandbox_id=sandbox.object_id,image_id=IMAGE)
        statefile.write_text(json.dumps(state,indent=2)+'\n');print('STARTED',sandbox.object_id,flush=True);return
    sandbox=modal.Sandbox.from_id(state['sandbox_id']);rc=sandbox.poll();print('POLL',sandbox.object_id,rc,flush=True)
    if rc is None:return
    (out/'stdout.log').write_text(sandbox.stdout.read());(out/'stderr.log').write_text(sandbox.stderr.read());state['exit_code']=rc;statefile.write_text(json.dumps(state,indent=2)+'\n')
    if rc:raise RuntimeError((out/'stderr.log').read_text()[-5000:])
    prefix=remote.removeprefix('/data/')+'/results/'
    receipt=json.loads(b''.join(volume.read_file(prefix+'complete.json')))
    assert receipt['request_sha256']==request_sha and receipt['worker_sha256']==sha(worker)
    for name,r in receipt['outputs'].items():
        payload=b''.join(volume.read_file(prefix+name));assert len(payload)==r['bytes'] and hashlib.sha256(payload).hexdigest()==r['sha256']
        (out/name).write_bytes(payload)
    (out/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n');state['status']='COMPLETED';statefile.write_text(json.dumps(state,indent=2)+'\n')
    verify_outputs(out);print('COLLECTED',out,flush=True)
    from pool_report import replay_and_report
    print(json.dumps(replay_and_report(out),indent=2))


if __name__=='__main__':main()
