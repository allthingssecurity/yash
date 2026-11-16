"""Plan or execute a read-only RNA-difference/normal-self audit on Modal."""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2]
VOLUME='research-lens-workspace'
PARENT='/data/results/TN25-112103/matched-normal-20261107'
IMAGE='im-hX1h6HknDMvasm9Z6q530H'
GENES={'F5','MUC5B','MUC3A','MUC6'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design',type=Path,default=Path('designs/design-3'))
    p.add_argument('--output',type=Path,default=Path('work/rna-normal-self-audit-20261112'))
    p.add_argument('--execute',action='store_true');p.add_argument('--collect',action='store_true')
    p.add_argument('--run-id',default='rna-normal-self-audit-20261112',help='Use a new name for a fresh remote rerun; completed measurements are never overwritten')
    a=p.parse_args();assert re.fullmatch('[a-z0-9-]+',a.run_id), 'Use a simple run identifier'
    design=(ROOT/a.design).resolve();out=(ROOT/a.output).resolve()
    pool=design/'inputs/sequence-review/candidate-pool.json'
    if not pool.exists():pool=design/'selection/selected.json'
    selected=json.loads(pool.read_text())
    chosen={r['window_id']:r for r in selected if r['source']=='SPLICE' and r['gene'] in GENES}
    assert len(chosen)==4
    source=ROOT/'data/shared-rna/splicing/reconstruction/patient-windows.jsonl.gz'
    verification=ROOT/'data/shared-rna/splicing/read-verification/verification.json'
    v=json.loads(verification.read_text());assert v['status']=='PASS' and sha(source)==v['reconstruction_sha256']
    windows=[]
    with gzip.open(source,'rt') as f:
        for line in f:
            w=json.loads(line)
            if w['window_id'] in chosen:
                r=chosen[w['window_id']];assert w['sequence']==r['sequence']
                w['region_id']=r['region_id'];windows.append(w)
    assert len(windows)==4
    ids={w['region_id'] for w in windows}
    additional=[]
    for r in selected:
        if r['source'] not in ['MUTATION','SPLICE','FUSION'] or r['region_id'] in ids:continue
        row={k:r[k] for k in ['gene','region_id','sequence','source']}
        if 'target_start_in_window' in r:
            row.update(target_start_in_window=int(r['target_start_in_window']),target_end_in_window=int(r['target_end_in_window']))
        else:
            assert r['context_id']=='CDH1_NRG1_INFRAME'
            fusion=json.loads((ROOT/'data/shared-rna/fusions/inframe-request.json').read_text())
            junction=fusion['junction_codon_0based']-r['window_start']
            assert r['sequence']==fusion['source_sequence'][r['window_start']:r['window_end']]
            row.update(target_start_in_window=junction,target_end_in_window=junction+1)
        additional.append(row)
    request=dict(design=str(design.relative_to(ROOT)),selected_sha256=sha(pool),
        patient_window_sha256=sha(source),existing_BAM_verification_sha256=sha(verification),windows=windows,additional_regions=additional,
        scope='Fresh matched-normal/tumour focal-base checks for four RNA-difference regions; broad reference/personal-self search for all selected altered-source regions; CTA/self excluded')
    body=(json.dumps(request,indent=2)+'\n').encode()
    request_sha=hashlib.sha256(body).hexdigest()
    worker=ROOT/'tools/rna_evidence/worker.py';bundle=[worker,ROOT/'designs/design-2/pipeline/normal_evidence.py',ROOT/'designs/design-2/pipeline/self_search.py']
    digest=hashlib.sha256(b''.join(p.name.encode()+p.read_bytes() for p in bundle)).hexdigest()
    remote=f'/data/results/TN25-112103/{a.run_id}/{request_sha[:12]}-{digest[:12]}'
    print(json.dumps(dict(remote=remote,design=str(design),normal_base_review_genes=sorted(GENES),altered_regions=len(windows)+len(additional),execute=a.execute,collect=a.collect),indent=2),flush=True)
    if not a.execute and not a.collect:return
    import modal
    out.mkdir(parents=True,exist_ok=True);statefile=out/'state.json';reqfile=out/'request.json'
    if reqfile.exists():assert reqfile.read_bytes()==body,'Use a new output for different inputs'
    else:reqfile.write_bytes(body)
    vol=modal.Volume.from_name(VOLUME)
    state=json.loads(statefile.read_text()) if statefile.exists() else dict(remote=remote,request_sha256=request_sha,bundle_sha256=digest)
    assert state['remote']==remote
    if 'sandbox_id' not in state:
        if a.collect:raise RuntimeError('Audit not started; --execute first')
        for source_path,key in [(reqfile,remote+'/request.json')]+[(p,remote+'/code/'+p.name) for p in bundle]:
            try:
                with vol.batch_upload() as upload:upload.put_file(str(source_path),key.removeprefix('/data/'))
            except FileExistsError:assert b''.join(vol.read_file(key.removeprefix('/data/')))==source_path.read_bytes()
        sandbox=modal.Sandbox.create('python','-u',remote+'/code/worker.py','--root',PARENT,'--request',remote+'/request.json','--output',remote+'/results',
            app=modal.App.lookup('research-yash-matched-normal'),image=modal.Image.from_id(IMAGE),cpu=2,memory=16384,timeout=7200,volumes={'/data':vol})
        state['sandbox_id']=sandbox.object_id;state['image_id']=IMAGE
        statefile.write_text(json.dumps(state,indent=2)+'\n')
        print('STARTED',sandbox.object_id,flush=True)
        return
    sandbox=modal.Sandbox.from_id(state['sandbox_id']);rc=sandbox.poll()
    print('POLL',sandbox.object_id,rc,flush=True)
    if rc is None:return
    (out/'stdout.log').write_text(sandbox.stdout.read());(out/'stderr.log').write_text(sandbox.stderr.read())
    state['exit_code']=rc;statefile.write_text(json.dumps(state,indent=2)+'\n')
    if rc:raise RuntimeError((out/'stderr.log').read_text()[-5000:])
    key=remote.removeprefix('/data/')+'/results/'
    receipt=json.loads(b''.join(vol.read_file(key+'complete.json')))
    assert receipt['status']=='COMPLETED' and receipt['request_sha256']==request_sha and receipt['worker_sha256']==sha(worker)
    for name,record in receipt['outputs'].items():
        content=b''.join(vol.read_file(key+name));assert len(content)==record['bytes'] and hashlib.sha256(content).hexdigest()==record['sha256']
        (out/name).write_bytes(content)
    (out/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n')
    state['status']='COMPLETED';statefile.write_text(json.dumps(state,indent=2)+'\n')
    print('COLLECTED',out,flush=True)


if __name__=='__main__':main()
