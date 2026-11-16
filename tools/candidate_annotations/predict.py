"""Fresh, context-keyed MHCflurry inference inside the pinned Modal runtime."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import io
import json
import math
import os
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--request', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--models', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    os.environ['OMP_NUM_THREADS'] = '4'
    os.environ['MKL_NUM_THREADS'] = '4'
    import numpy as np
    import torch
    from mhcflurry import Class1PresentationPredictor
    torch.set_num_threads(4)
    torch.set_float32_matmul_precision('highest')
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    np.random.seed(0)
    torch.manual_seed(0)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(0)
    req = json.loads(a.request.read_text())
    expected = json.loads(a.manifest.read_text())
    observed = [dict(path=str(p.relative_to(a.models)), bytes=p.stat().st_size, sha256=sha(p))
                for p in sorted(a.models.rglob('*')) if p.is_file()]
    if observed != expected:
        raise ValueError('Model bundle differs from the pinned manifest')
    versions = {k:importlib.metadata.version(k) for k in ['mhcflurry','numpy','torch']}
    if versions['mhcflurry'] != '2.3.14':
        raise ValueError('The annotation stage pins MHCflurry 2.3.14')
    model = Class1PresentationPredictor.load(str(a.models))
    if not model.supports_presentation_prediction:
        raise ValueError('Presentation components are unavailable')
    lengths = {k:int(v) for k,v in model.processing_predictor_with_flanks.sequence_lengths.items()}
    if lengths != dict(peptide=15,n_flank=5,c_flank=5):
        raise ValueError('Unexpected processing-model context lengths')
    started = time.perf_counter()
    rows = []
    for mode in ['without_flanks','with_flanks']:
        queries = [q for q in req['queries'] if q['mode'] == mode]
        if not queries:
            continue
        arguments = dict(peptides=[q['peptide'] for q in queries], include_affinity_percentile=True, verbose=0)
        if mode == 'with_flanks':
            arguments.update(n_flanks=[q['n_flank'] for q in queries],c_flanks=[q['c_flank'] for q in queries])
        for allele in req['alleles']:
            frame = model.predict(alleles=[allele], **arguments)
            if frame['peptide'].tolist() != arguments['peptides']:
                raise ValueError('Prediction rows changed query order')
            for query,value in zip(queries,frame.to_dict('records')):
                row = dict(query_id=query['query_id'],peptide=query['peptide'],n_flank=query['n_flank'],
                    c_flank=query['c_flank'],mode=mode,allele=allele,
                    affinity_nm=float(value['affinity']),affinity_percentile=float(value['affinity_percentile']),
                    processing_score=float(value['processing_score']),presentation_score=float(value['presentation_score']))
                if not all(math.isfinite(row[k]) for k in ['affinity_nm','affinity_percentile','processing_score','presentation_score']):
                    raise ValueError('Non-finite model prediction')
                if row['affinity_nm'] <= 0 or not all(0 <= row[k] <= 1 for k in ['processing_score','presentation_score']):
                    raise ValueError('Prediction outside model bounds')
                rows.append(row)
            print('PREDICTED',mode,allele,len(queries),flush=True)
    a.output.mkdir(parents=True,exist_ok=False)
    stream = io.StringIO()
    writer = csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    writer.writeheader()
    writer.writerows(sorted(rows,key=lambda r:(r['query_id'],r['allele'])))
    with (a.output/'predictions.tsv.gz').open('wb') as f:
        with gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as z:
            z.write(stream.getvalue().encode())
    meta = dict(status='COMPLETED',new_model_inference=True,model_files_verified=True,
        model_manifest_sha256=sha(a.manifest),model_files=len(observed),request_sha256=sha(a.request),
        prediction_sha256=sha(a.output/'predictions.tsv.gz'),worker_sha256=sha(Path(__file__)),
        versions=versions,GPU=torch.cuda.get_device_name() if torch.cuda.is_available() else None,
        elapsed_seconds=time.perf_counter()-started,queries=len(req['queries']),prediction_rows=len(rows),
        context_is_part_of_query_identity=True,uses_RNA_count_or_expression_as_model_input=False,
        models_redistributed=False)
    (a.output/'inference.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(meta),flush=True)


if __name__ == '__main__':
    main()
