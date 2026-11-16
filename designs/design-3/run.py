"""Plan, launch or collect the final 30-region construct build on Modal."""
from pathlib import Path
import argparse
import asyncio
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
REMOTE='/data/results/TN25-112103/design-3-final-20261116/attempt-01'
LOCAL=REPO/'work/design-3-final-20261116/attempt-01'
VOLUME='research-lens-workspace'
PHASES=['prepare','predict','solve','predict-selected','verify','encode']
FOLDERS=dict(prepare='prepare',predict='all-inference',solve='protein',**{'predict-selected':'selected-inference'},verify='verify',encode='mrna')
IMAGES=dict(prepare='im-hX1h6HknDMvasm9Z6q530H',predict='im-DBMghmAQfB3Oz7McdNZeic',
    solve='im-bgDqBGTPmD1LPQS5T4tfm8',**{'predict-selected':'im-DBMghmAQfB3Oz7McdNZeic'},
    verify='im-hX1h6HknDMvasm9Z6q530H',encode='im-yWdXWgFXTigt4nutlA984S')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path):return json.loads(Path(path).read_text())
def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp');temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def inputs():
    mapping={f'code/{p.name}':p for p in sorted((HERE/'pipeline').glob('*.py'))}
    mapping['policy.json']=HERE/'policy.json'
    mapping.update({'integrate/'+p.name:p for p in sorted((HERE/'selection').iterdir()) if p.is_file()})
    mapping['case.json']=HERE/'inputs/case.json'
    for name in ['__init__.py','hydrophobicity.py','io.py','paths.py']:
        mapping['frozen/src/yash_vaccine/'+name]=REPO/'src/yash_vaccine'/name
    for path in ['config/linkers.json','designs/initial/construct/architecture.json',
                 'designs/initial/construct/predict.py','designs/initial/construct/optimize.py']:
        mapping['frozen/'+path]=REPO/path
    for name in ['run.py','optimizer.py','genbank.py','backbone.json','annotation-sources.json','human-codon-usage.json','requirements.txt']:
        path='designs/initial/construct/mrna/'+name;mapping['frozen/'+path]=REPO/path
    for name in ['construct-inference.json','construct-model-files.json']:
        mapping['frozen/data/checkpoints/provenance/'+name]=REPO/'designs/initial/inputs/checkpoints/provenance'/name
    return mapping


def collect(volume,stage,state):
    record=state['stages'][stage]
    import modal
    sandbox=modal.Sandbox.from_id(record['sandbox_id']);rc=sandbox.poll();record['last_poll_exit_code']=rc
    print('POLL',stage,record['sandbox_id'],rc,flush=True)
    if rc is None:return False
    if rc!=0:
        record['status']='FAILED'
        async def stderr():return await asyncio.wait_for(sandbox.stderr.read.aio(),timeout=5)
        try:
            log=asyncio.run(stderr());(LOCAL/(stage+'.stderr.log')).write_text(log);print(log[-4000:],flush=True)
        except Exception as error:record['log_observation_error']=type(error).__name__
        save(LOCAL/'state.json',state)
        raise RuntimeError(stage+' has a confirmed nonzero exit; inspect the retained worker before retrying')
    folder=FOLDERS[stage];prefix=REMOTE.removeprefix('/data/')+'/construct/'+folder+'/'
    receipt=json.loads(b''.join(volume.read_file(prefix+'complete.json')))
    assert receipt['status']=='COMPLETED'
    target=HERE/'construct'/folder;target.mkdir(parents=True,exist_ok=True)
    for name,r in receipt['outputs'].items():
        if r['bytes']>40_000_000:continue
        p=target/name
        if p.exists() and p.stat().st_size==r['bytes'] and sha(p)==r['sha256']:continue
        body=b''.join(volume.read_file(prefix+name))
        assert len(body)==r['bytes'] and hashlib.sha256(body).hexdigest()==r['sha256'],name
        p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(body)
    save(target/'complete.json',receipt)
    if stage=='prepare':
        body=b''.join(volume.read_file(REMOTE.removeprefix('/data/')+'/construct/architecture.json'))
        architecture=json.loads(body);expected=load(REPO/'designs/initial/construct/architecture.json');expected['name']=load(HERE/'policy.json')['name'];assert architecture==expected
        (HERE/'construct/architecture.json').write_bytes(body)
    if stage=='verify':
        (HERE/'construct/final-check.json').write_bytes((target/'final-check.json').read_bytes())
    if stage=='encode':
        prefix=REMOTE.removeprefix('/data/')+'/rebuild/construct/mrna/'
        for name in ['run.py','optimizer.py','genbank.py','backbone.json','annotation-sources.json','human-codon-usage.json','requirements.txt']:
            (target/name).write_bytes(b''.join(volume.read_file(prefix+name)))
    record.update(status='COMPLETED',receipt_sha256=sha(target/'complete.json'))
    save(LOCAL/'state.json',state);print('COLLECTED',stage,flush=True)
    return True


def main():
    p=argparse.ArgumentParser(description=__doc__);mode=p.add_mutually_exclusive_group()
    mode.add_argument('--execute',action='store_true');mode.add_argument('--collect',action='store_true')
    mode.add_argument('--replay',action='store_true');p.add_argument('--output',type=Path);p.add_argument('--through',choices=PHASES)
    p.add_argument('--worker-revision',action='store_true',help='Preserve the first deployed worker and use a separately frozen worker bundle for later stages')
    a=p.parse_args()
    subprocess.run([sys.executable,str(HERE/'selection.py'),'--verify'],check=True)
    if a.replay:
        command=[sys.executable,str(HERE/'replay.py')]
        if a.output:command+=['--output',str(a.output)]
        subprocess.run(command,check=True);return
    files=inputs();manifest={path:dict(repository_path=str(local.relative_to(REPO)),bytes=local.stat().st_size,sha256=sha(local)) for path,local in files.items()}
    digest=hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest()
    print(json.dumps(dict(execute=a.execute,collect=a.collect,regions=30,stages=PHASES,remote=REMOTE,source_bundle_sha256=digest,
        scope='Fresh joint assembly and RNA encoding of the agreed final 30-region selection; no new variant calling, lead promotion, publication or website change'),indent=2),flush=True)
    if not (a.execute or a.collect):return
    import modal
    LOCAL.mkdir(parents=True,exist_ok=True);statepath=LOCAL/'state.json'
    state=load(statepath) if statepath.exists() else dict(root=REMOTE,source_bundle_sha256=digest,stages={})
    volume=modal.Volume.from_name(VOLUME)
    if state['source_bundle_sha256']!=digest:
        assert a.execute and a.worker_revision,'Worker changes require an explicit preserved source revision'
        previous=load(HERE/'execution/source-v1/inputs.manifest.json')
        assert sha(HERE/'execution/source-v1/cassette.py')==previous['code/cassette.py']['sha256']
        assert set(previous)==set(manifest)
        assert all(previous[k]==v for k,v in manifest.items() if not k.startswith('code/')),'Scientific inputs cannot change inside a build'
        prefix='code/revisions/'+digest[:16]+'/'
        revised={prefix+Path(k).name if k.startswith('code/') else k:v for k,v in manifest.items()}
        native_manifest=REMOTE+'/worker-manifests/'+digest+'.json'
        save(LOCAL/('inputs-'+digest+'.json'),revised)
        with volume.batch_upload() as upload:
            for k,path in files.items():
                if k.startswith('code/'):upload.put_file(str(path),(REMOTE+'/'+prefix+path.name).removeprefix('/data/'))
            upload.put_file(str(LOCAL/('inputs-'+digest+'.json')),native_manifest.removeprefix('/data/'))
        original=state['source_bundle_sha256'];state.setdefault('previous_source_bundles',[]).append(original)
        for record in state['stages'].values():record.setdefault('source_bundle_sha256',original)
        state.update(source_bundle_sha256=digest,worker=REMOTE+'/'+prefix+'cassette.py',worker_manifest=native_manifest)
        save(HERE/'execution/inputs.manifest.json',revised);save(statepath,state)
    if a.execute and not state.get('uploaded'):
        save(LOCAL/'inputs.manifest.json',manifest)
        with volume.batch_upload() as upload:
            for destination,path in files.items():upload.put_file(str(path),(REMOTE+'/'+destination).removeprefix('/data/'))
            upload.put_file(str(LOCAL/'inputs.manifest.json'),(REMOTE+'/inputs.manifest.json').removeprefix('/data/'))
        state['uploaded']=True;save(statepath,state)
        save(HERE/'execution/inputs.manifest.json',manifest)
    for stage in PHASES:
        previous=state['stages'].get(stage)
        if previous:
            if previous['status']=='FAILED':raise RuntimeError('Inspect confirmed failed stage before creating a new attempt: '+stage)
            if previous['status']!='COMPLETED' and not collect(volume,stage,state):save(statepath,state);return
        elif not a.execute:return
        else:
            gpu=stage in ['predict','predict-selected']
            with modal.enable_output():
                command=['python','-u',state.get('worker',REMOTE+'/code/cassette.py'),REMOTE,stage]
                if state.get('worker_manifest'):command.append(state['worker_manifest'])
                sb=modal.Sandbox.create(*command,
                    app=modal.App.lookup('research-yash-matched-normal'),image=modal.Image.from_id(IMAGES[stage]),
                    gpu='A10' if gpu else None,cpu=4 if gpu else 2,memory=32768 if gpu else 16384,
                    timeout=21600,volumes={'/data':volume})
            state['stages'][stage]=dict(sandbox_id=sb.object_id,status='RUNNING',image_id=IMAGES[stage],gpu='A10' if gpu else None,source_bundle_sha256=digest)
            save(statepath,state);print('STARTED',stage,sb.object_id,flush=True);return
        if a.through==stage:return
    state['status']='COMPLETED';save(statepath,state);save(HERE/'execution/modal.json',state)
    print('DESIGN3_FINAL_BUILD_COMPLETED',flush=True)


if __name__=='__main__':main()
