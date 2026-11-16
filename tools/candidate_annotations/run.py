"""Plan, execute or offline-verify candidate annotations; no selection mutation."""
from pathlib import Path
import argparse
import datetime
import json

from analysis import ROOT, DEFAULT, request, json_bytes, sha, finish

MODELS = '/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=DEFAULT)
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--execute', action='store_true')
    mode.add_argument('--verify', action='store_true')
    a = p.parse_args()
    out = a.output.resolve()
    if a.verify:
        print(json.dumps(finish(out,verify=True),indent=2))
        return
    req = request()
    print(json.dumps(dict(regions=len(req['regions']),current_plan_regions=30,
        queries=len(req['queries']),peptide_HLA_predictions=len(req['queries'])*6,
        execute=a.execute,selection_changed=False),indent=2),flush=True)
    if not a.execute:
        return
    if (out/'request.json').exists():
        raise ValueError('Choose a fresh output directory; existing inference is immutable')
    import modal
    out.mkdir(parents=True,exist_ok=True)
    (out/'request.json').write_bytes(json_bytes(req))
    worker = Path(__file__).with_name('predict.py')
    manifest = ROOT/'designs/initial/inputs/checkpoints/provenance/construct-model-files.json'
    remote = '/tmp/yash-candidate-annotations'
    volume = modal.Volume.from_name('research-lens-workspace')
    state = dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        request_sha256=sha(out/'request.json'),worker_sha256=sha(worker),terminated=False)
    sandbox = None
    try:
        sandbox = modal.Sandbox.create(app=modal.App.lookup('research-yash-vaccine',create_if_missing=True),
            image=modal.Image.from_id('im-DBMghmAQfB3Oz7McdNZeic'),gpu='A10',cpu=4,memory=16384,
            timeout=1800,volumes={'/data':volume})
        state['sandbox_id'] = sandbox.object_id
        (out/'execution.json').write_bytes(json_bytes(state))
        for path,name in [(out/'request.json','request.json'),(worker,'predict.py'),(manifest,'model-manifest.json')]:
            sandbox.filesystem.write_bytes(path.read_bytes(),remote+'/'+name)
        proc = sandbox.exec('python','-u',remote+'/predict.py','--request',remote+'/request.json',
            '--manifest',remote+'/model-manifest.json','--models',MODELS,'--output',remote+'/output',timeout=1500)
        for line in proc.stdout:
            print(line,end='',flush=True)
        error = proc.stderr.read()
        state['exit_code'] = proc.wait()
        if error:
            (out/'stderr.log').write_text(error)
        if state['exit_code']:
            raise RuntimeError('Prediction failed: '+error[-5000:])
        for name in ['predictions.tsv.gz','inference.json']:
            (out/name).write_bytes(sandbox.filesystem.read_bytes(remote+'/output/'+name))
        inference = json.loads((out/'inference.json').read_text())
        if inference['worker_sha256'] != sha(worker):
            raise ValueError('Remote worker identity mismatch')
        print(json.dumps(finish(out),indent=2),flush=True)
    finally:
        if sandbox:
            sandbox.terminate(wait=True)
            state['terminated'] = True
        state['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        (out/'execution.json').write_bytes(json_bytes(state))


if __name__ == '__main__':
    main()
